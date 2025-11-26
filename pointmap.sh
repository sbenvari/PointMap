#!/usr/bin/env bash
set -euo pipefail

###############################################
# gene_pipeline.sh
# Extract gene → annotate → extract proteins → add reference → align → call mutations
# USAGE:
#   ./gene_pipeline.sh <GENOMES_DIR> <REF_GENE.fasta> <GENE_NAME> <OUTPUT_DIR>
###############################################

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <GENOMES_DIR> <REF_GENE.fasta> <GENE_NAME> <OUTPUT_DIR>"
    exit 1
fi

GENOMES_DIR="$1"
REF_GENE="$2"
GENE_NAME="$3"
OUT="$4"

mkdir -p "$OUT"
SEQ="$OUT/01_sequences"
PROKKA="$OUT/02_prokka"
PROT="$OUT/03_proteins"
ALIGN="$OUT/04_alignment"
RESULTS="$OUT/05_results"

mkdir -p "$SEQ" "$PROKKA" "$PROT" "$ALIGN" "$RESULTS"

shopt -s nullglob

###############################################
# STEP 0: Convert REFERENCE gene to a REF.faa
###############################################

echo "=== Step 0: Processing reference gene ==="

mkdir -p "$OUT/reference_tmp"

prokka \
  --quiet \
  --locustag REF \
  --outdir "$OUT/reference_tmp" \
  "$REF_GENE"

REF_FAA=$(find "$OUT/reference_tmp" -name "*.faa" | head -n 1)

if [[ ! -f "$REF_FAA" ]]; then
    echo "ERROR: reference protein .faa not produced by Prokka"
    exit 1
fi

# Clean + rename reference protein
awk -v id="reference" 'NR==1 {print ">"id; next} {print}' "$REF_FAA" \
    > "$PROT/reference.faa"

echo "Reference protein created: $PROT/reference.faa"


###############################################
# STEP 1: Extract gene from sample genomes
###############################################

echo "=== Step 1: Extracting $GENE_NAME from genomes ==="

for genome in "$GENOMES_DIR"/*.fa*; do
    sample=$(basename "$genome")
    sample="${sample%%.*}"

    blastn -query "$REF_GENE" -subject "$genome" -out "${sample}.blast" -outfmt 6

    if [[ ! -s ${sample}.blast ]]; then
        echo "WARNING: No BLAST hit for $sample"
        continue
    fi

    awk '{ if ($9 < $10) print $2 "\t" $9-1 "\t" $10;
           else           print $2 "\t" $10-1 "\t" $9 }' \
        "${sample}.blast" > "${sample}.bed"

    bedtools getfasta -fi "$genome" -bed "${sample}.bed" \
        -fo "$SEQ/${sample}_${GENE_NAME}.fasta"

    rm -f "${sample}.blast" "${sample}.bed"
done


###############################################
# STEP 2: Clean FASTA headers
###############################################

echo "=== Step 2: Cleaning FASTA headers ==="

for f in "$SEQ"/*.fasta; do
    awk '/^>/ {sub(/_length_.*/, "", $0); print} !/^>/ {print}' "$f" \
        > tmp && mv tmp "$f"
done


###############################################
# STEP 3: Prokka annotation of extracted sequences
###############################################

echo "=== Step 3: Running Prokka ==="

for f in "$SEQ"/*.fasta; do
    sample=$(basename "$f" .fasta)
    prokka --quiet --outdir "$PROKKA/$sample" "$f"
done


###############################################
# STEP 4: Extract protein sequences
###############################################

echo "=== Step 4: Extracting proteins ==="

for d in "$PROKKA"/*/; do
    sample=$(basename "$d")
    faa=$(find "$d" -name "*.faa" | head -n 1)

    if [[ -f "$faa" ]]; then
        awk -v id="$sample" 'NR==1 {print ">" id; next} {print}' "$faa" \
            > "$PROT/${sample}.faa"
    fi
done


###############################################
# STEP 5: Merge reference + all protein sequences
###############################################

echo "=== Step 5: Building multi-protein file ==="

cat "$PROT/reference.faa" "$PROT"/*.faa > "$ALIGN/all.faa"


###############################################
# STEP 6: MAFFT alignment
###############################################

echo "=== Step 6: Running MAFFT ==="
mafft --auto "$ALIGN/all.faa" > "$ALIGN/aligned.faa"


###############################################
# STEP 7: Python mutation caller
###############################################

echo "=== Step 7: Mutation calling ==="

python3 - <<EOF
from Bio import AlignIO

alignment = AlignIO.read("$ALIGN/aligned.faa", "fasta")

ref = alignment[0]
ref_seq = ref.seq

with open("$RESULTS/${GENE_NAME}_mutations.txt", "w") as out:
    for record in alignment[1:]:
        muts = []
        for i, (r,s) in enumerate(zip(ref_seq, record.seq), start=1):
            if r != s and r != "-" and s != "-":
                muts.append(f"{r}{i}{s}")
        if muts:
            out.write(f"{record.id}\t{','.join(muts)}\n")
        else:
            out.write(f"{record.id}\tNo mutations\n")
EOF

echo "=== Pipeline completed ==="
echo "Results stored in $OUT"
