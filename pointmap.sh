#!/usr/bin/env bash
set -euo pipefail

# USAGE:
#   ./pointmap.sh <REFERENCE_GENOME> <GENOMES_DIR> <GENE_NAME> <OUTDIR>

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <REFERENCE_GENOME> <GENOMES_DIR> <GENE_NAME> <OUTDIR>"
    exit 1
fi

# Resolve absolute paths
REF_GENOME=$(realpath "$1")
SAMPLE_DIR=$(realpath "$2")
GENE="$3"
OUT=$(realpath "$4")

mkdir -p "$OUT"
cd "$OUT"

echo "=== POINTMAP PIPELINE STARTED ==="

############################################
# GLOBAL CLEANUP OF .fai FILES
############################################
echo "=== CLEANUP: Removing .fai index files ==="
find "$SAMPLE_DIR" -maxdepth 1 -name "*.fai" -delete 2>/dev/null || true
find "$OUT" -maxdepth 4 -name "*.fai" -delete 2>/dev/null || true

############################################
# STEP 1: Annotate REFERENCE genome
############################################
echo "=== STEP 1: Running Prokka on reference genome ==="
prokka --force --prefix ref --outdir ref_prokka "$REF_GENOME"

############################################
# STEP 2: Extract gene coordinates from GFF
############################################
echo "=== STEP 2: Locating $GENE in reference annotation ==="

grep -i "Name=${GENE}" ref_prokka/*.gff > ref_${GENE}.gff

if [[ ! -s ref_${GENE}.gff ]]; then
    echo "ERROR: Gene '$GENE' not found in reference annotation."
    exit 1
fi

awk '{print $1"\t"$4-1"\t"$5}' ref_${GENE}.gff > ref_${GENE}.bed

############################################
# STEP 3: Extract reference gene FASTA
############################################
echo "=== STEP 3: Extracting reference $GENE sequence ==="

mkdir -p sequences_ref

bedtools getfasta -fi "$REF_GENOME" -bed ref_${GENE}.bed \
    -fo sequences_ref/${GENE}_reference.fasta

sed -i -e '$a\' sequences_ref/${GENE}_reference.fasta

############################################
# STEP 4: Extract gene from sample genomes
############################################
echo "=== STEP 4: Extracting $GENE from sample genomes ==="

mkdir -p sequences_samples

for genome in "$SAMPLE_DIR"/*.fa "$SAMPLE_DIR"/*.fna "$SAMPLE_DIR"/*.fasta; do
    [[ -e "$genome" ]] || continue
    [[ "$genome" == *.fai ]] && continue

    [[ -e "${genome}.fai" ]] && rm -f "${genome}.fai"

    sample=$(basename "$genome")
    sample="${sample%%.*}"

    blastn -query sequences_ref/${GENE}_reference.fasta \
           -subject "$genome" -out ${sample}.blast -outfmt 6 || true

    if [[ ! -s ${sample}.blast ]]; then
        echo "WARNING: No BLAST hit for $sample"
        continue
    fi

    awk '{if ($9 < $10) print $2"\t"$9-1"\t"$10;
          else print $2"\t"$10-1"\t"$9}' \
        ${sample}.blast > ${sample}.bed

    bedtools getfasta -fi "$genome" -bed ${sample}.bed \
        -fo sequences_samples/${sample}_${GENE}.fasta

    rm -f ${sample}.blast ${sample}.bed
    rm -f "${genome}.fai" 2>/dev/null || true
done

############################################
# STEP 5: Rename FASTA headers to EXACT sample names
############################################
echo "=== STEP 5: Renaming FASTA headers to sample names ==="

# Reference header
awk 'BEGIN {print ">reference"} !/^>/ {print}' \
    sequences_ref/${GENE}_reference.fasta > tmp && \
mv tmp sequences_ref/${GENE}_reference.fasta

# Sample headers
for f in sequences_samples/*.fasta; do
    sample=$(basename "$f" .fasta)

    awk -v id="$sample" 'BEGIN {print ">" id} !/^>/ {print}' \
        "$f" > tmp && mv tmp "$f"
done

############################################
# STEP 6: Prokka → proteins
############################################
echo "=== STEP 6: Running Prokka ==="

mkdir -p prokka_out protein

# Reference
prokka --quiet --force --prefix reference --outdir prokka_out/reference \
       sequences_ref/${GENE}_reference.fasta

cp prokka_out/reference/reference.faa protein/reference.faa

# Samples
for f in sequences_samples/*.fasta; do
    sample=$(basename "$f" .fasta)
    outdir="prokka_out/$sample"

    prokka --quiet --force --prefix "$sample" --outdir "$outdir" "$f"

    faa="$outdir/${sample}.faa"
    if [[ -f "$faa" ]]; then
        cp "$faa" protein/${sample}.faa
    else
        echo "WARNING: No protein file for $sample"
    fi
done

############################################
# STEP 6.5: Normalize protein FASTA headers
############################################
echo "=== STEP 6.5: Renaming protein headers to sample names ==="

for faa in protein/*.faa; do
    sample=$(basename "$faa" .faa)

    awk -v id="$sample" '
        /^>/ {print ">" id; next}
        {print}
    ' "$faa" > tmp && mv tmp "$faa"
done

############################################
# STEP 7: Concatenate proteins
############################################
echo "=== STEP 7: Concatenating FAA files ==="
cat protein/*.faa > multi_sequences.faa

############################################
# STEP 8: MAFFT alignment
############################################
echo "=== STEP 8: Running MAFFT ==="
mafft --auto multi_sequences.faa > aligned_sequences.faa

############################################
# STEP 9: Mutation calling (direct sample names)
############################################
############################################
# STEP 9: Mutation calling (reference excluded)
############################################
echo "=== STEP 9: Calling mutations ==="

python3 - <<EOF
from Bio import AlignIO

alignment = AlignIO.read("aligned_sequences.faa", "fasta")

# Identify the reference record by its header (case-insensitive)
ref_id_candidates = ["reference", "ref", "REF"]
ref_index = None
for i, rec in enumerate(alignment):
    if rec.id in ref_id_candidates:
        ref_index = i
        break

if ref_index is None:
    raise ValueError("ERROR: Reference sequence not found in alignment.")

ref = alignment[ref_index]
ref_seq = ref.seq

with open("mutations_${GENE}.txt", "w") as out:
    out.write("Sample\tMutations\n")

    for i, record in enumerate(alignment):
        if i == ref_index:
            continue  # <-- SKIP REFERENCE

        sid = record.id
        muts = []

        for pos, (r, s) in enumerate(zip(ref_seq, record.seq), start=1):
            if r != s and r != '-' and s != '-':
                muts.append(f"{r}{pos}{s}")

        out.write(f"{sid}\t{','.join(muts) if muts else 'No mutations'}\n")
EOF
