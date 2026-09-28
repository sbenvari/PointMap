# PointMap — Nextflow version

A Nextflow DSL2 reimplementation of the original PointMap Bash workflow for extracting a named gene from a reference and a set of sample genomes, translating/annotating the extracted sequences, aligning proteins, and reporting amino-acid substitutions relative to the reference.

## Workflow

1. Annotate the reference genome with Prokka.
2. Locate the requested gene in the reference GFF.
3. Extract the reference gene with BEDTools.
4. BLAST the reference gene against each sample genome in parallel.
5. Extract the corresponding sequence from each sample genome.
6. Annotate extracted gene sequences with Prokka and collect proteins.
7. Concatenate proteins.
8. Align proteins with MAFFT.
9. Call amino-acid substitutions relative to the reference with Biopython.

## Requirements

- Nextflow
- Java (required by Nextflow)
- Conda/Mamba if using the supplied `conda` profile

The bioinformatics dependencies are listed in `environment.yml`.

## Run

```bash
nextflow run main.nf \
  --ref reference.fasta \
  --genomes genomes \
  --gene gyrA \
  --outdir results \
  -profile conda
```

`--genomes` should point to a directory containing `.fa`, `.fna`, or `.fasta` files.

If your cluster uses SLURM:

```bash
nextflow run main.nf \
  --ref reference.fasta \
  --genomes genomes \
  --gene gyrA \
  --outdir results \
  -profile slurm
```

Adjust the resource requests in `nextflow.config` to match your system.

## Main outputs

- `results/reference_annotation/` — Prokka reference annotation
- `results/reference_gene/` — reference gene coordinates and FASTA
- `results/sequences_samples/` — extracted sample gene sequences
- `results/protein/` — normalised protein FASTA files
- `results/multi_sequences.faa` — concatenated proteins
- `results/aligned_sequences.faa` — MAFFT alignment
- `results/mutations_<GENE>.txt` — amino-acid substitutions by sample

## Important note on BLAST hits

This version intentionally preserves the logic of the original Bash workflow: all BLAST hits are converted to intervals and extracted. If a genome contains multiple hits, the current header-normalisation step collapses those extracted sequence records under one sample header, mirroring the original script's behaviour.

For a production version, a useful next improvement would be to define an explicit hit-selection rule (for example, best bit score / identity / coverage) and extract one orthologous hit per sample.

## Why Nextflow?

Compared with the original Bash implementation, this version makes each sample-level extraction an independent task, which allows Nextflow to parallelise work and manage intermediate files in isolated task directories. It also makes the workflow easier to resume, move to HPC, and package reproducibly with Conda.
