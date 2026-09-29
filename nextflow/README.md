# Nextflow implementation

The Nextflow version is intended for larger datasets, repeated analyses, or HPC use. Compared with the original Bash workflow, it can run sample-level tasks in parallel, resume completed work with `-resume`, isolate intermediate files, and automatically create a reproducible Conda environment from `environment.yml`.

For a small number of genomes, the Bash version may be simpler. For larger genome collections or cluster execution, the Nextflow version is more scalable and easier to reproduce.

## Requirements

- Conda
- Java 17 or later
- Nextflow

Create and activate a small environment for Nextflow and Java:

```bash
conda create -n nextflow -c conda-forge -c bioconda openjdk=17 nextflow
conda activate nextflow
```

Check the installation:

```bash
java -version
nextflow -version
```

The bioinformatics dependencies used by PointMap are defined in `environment.yml`. Nextflow creates and caches this analysis environment automatically during the first run, so it does not need to be installed manually.

## Usage

Run the workflow from the `nextflow/` directory:

```bash
nextflow run main.nf \
    --ref /path/to/reference.fasta \
    --genomes /path/to/genomes \
    --gene rpoB \
    --outdir /path/to/results
```

Example:

```bash
nextflow run main.nf \
    --ref ../example_ref/haemo_reference.fna \
    --genomes ../example_genome \
    --gene rpoB \
    --outdir ../rpoB-results
```

`--genomes` should point to a directory containing `.fa`, `.fna`, or `.fasta` genome files.

If a run is interrupted or a downstream step fails, rerun the same command with `-resume`:

```bash
nextflow run main.nf \
    --ref ../example_ref/haemo_reference.fna \
    --genomes ../example_genome \
    --gene rpoB \
    --outdir ../rpoB-results \
    -resume
```

Nextflow will reuse compatible completed tasks rather than repeating the full analysis.

The main final outputs are:

```text
multi_sequences.faa
aligned_sequences.faa
mutations_<GENE>.txt
```

Intermediate outputs, including extracted gene sequences and protein files, are also written to the selected output directory.
