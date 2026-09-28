#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * PointMap — Nextflow reimplementation of the original Bash workflow.
 *
 * Required parameters:
 *   --ref       Reference genome FASTA
 *   --genomes   Directory containing sample genomes (.fa/.fna/.fasta)
 *   --gene      Gene name as it appears in the Prokka GFF Name= attribute
 *
 * Optional:
 *   --outdir    Output directory (default: results)
 */

params.ref     = null
params.genomes = null
params.gene    = null
params.outdir  = 'results'

process ANNOTATE_REFERENCE {
    tag 'reference'
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}/reference_annotation", mode: 'copy'

    input:
    path ref_genome

    output:
    path 'ref_prokka/ref.gff', emit: gff

    script:
    """
    prokka \
        --force \
        --prefix ref \
        --outdir ref_prokka \
        ${ref_genome}
    """
}

process FIND_GENE {
    tag "${gene}"
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}/reference_gene", mode: 'copy'

    input:
    path gff
    val gene

    output:
    path "ref_${gene}.gff", emit: gene_gff
    path "ref_${gene}.bed", emit: bed

    script:
    """
    grep -i -F "Name=${gene}" ${gff} > ref_${gene}.gff || true

    if [[ ! -s ref_${gene}.gff ]]; then
        echo "ERROR: Gene '${gene}' not found in reference annotation." >&2
        exit 1
    fi

    awk 'BEGIN{OFS="\\t"} {print \$1, \$4-1, \$5}' \
        ref_${gene}.gff > ref_${gene}.bed
    """
}

process EXTRACT_REFERENCE_GENE {
    tag "${gene}"
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}/reference_gene", mode: 'copy'

    input:
    path ref_genome
    path bed
    val gene

    output:
    path "${gene}_reference.fasta", emit: fasta

    script:
    """
    bedtools getfasta \
        -fi ${ref_genome} \
        -bed ${bed} \
        -fo raw_reference.fasta

    # Preserve the original PointMap behaviour: one reference header.
    awk 'BEGIN {print ">reference"} !/^>/ {print}' \
        raw_reference.fasta > ${gene}_reference.fasta
    """
}

process EXTRACT_SAMPLE_GENE {
    tag "${sample}"
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}/sequences_samples", mode: 'copy'

    input:
    tuple val(sample), path(genome)
    path ref_gene
    val gene

    output:
    tuple val(sample), path("${sample}_${gene}.fasta"), emit: fasta, optional: true

    script:
    """
    blastn \
        -query ${ref_gene} \
        -subject ${genome} \
        -out ${sample}.blast \
        -outfmt 6 || true

    if [[ ! -s ${sample}.blast ]]; then
        echo "WARNING: No BLAST hit for ${sample}" >&2
        exit 0
    fi

    awk 'BEGIN{OFS="\\t"} {
        if (\$9 < \$10) print \$2, \$9-1, \$10;
        else             print \$2, \$10-1, \$9
    }' ${sample}.blast > ${sample}.bed

    bedtools getfasta \
        -fi ${genome} \
        -bed ${sample}.bed \
        -fo raw_${sample}.fasta

    # Preserve the original script's sample-header behaviour.
    awk -v id="${sample}_${gene}" \
        'BEGIN {print ">" id} !/^>/ {print}' \
        raw_${sample}.fasta > ${sample}_${gene}.fasta
    """
}

process ANNOTATE_GENE {
    tag "${sample}"
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}/protein", mode: 'copy'

    input:
    tuple val(sample), path(gene_fasta)

    output:
    tuple val(sample), path("${sample}.faa"), emit: protein

    script:
    """
    prokka \
        --quiet \
        --force \
        --prefix ${sample} \
        --outdir prokka_out \
        ${gene_fasta}

    faa="prokka_out/${sample}.faa"

    if [[ -s "\$faa" ]]; then
        awk -v id="${sample}" '
            /^>/ {print ">" id; next}
            {print}
        ' "\$faa" > ${sample}.faa
    else
        echo "WARNING: No protein sequence predicted for ${sample}" >&2
        touch ${sample}.faa
    fi
    """
}

process CONCAT_PROTEINS {
    tag 'all proteins'
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path faa_files

    output:
    path 'multi_sequences.faa', emit: multi

    script:
    def files = faa_files.join(' ')
    """
    cat ${files} > multi_sequences.faa
    """
}

process ALIGN_PROTEINS {
    tag 'MAFFT'
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path multi_faa

    output:
    path 'aligned_sequences.faa', emit: aligned

    script:
    """
    mafft --auto ${multi_faa} > aligned_sequences.faa
    """
}

process CALL_MUTATIONS {
    tag "${gene}"
    conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path alignment
    val gene

    output:
    path "mutations_${gene}.txt", emit: mutations

    script:
    """
    python3 ${projectDir}/bin/call_mutations.py \
        --alignment ${alignment} \
        --output mutations_${gene}.txt
    """
}

workflow {
    if (!params.ref) {
        error "Missing required parameter: --ref"
    }
    if (!params.genomes) {
        error "Missing required parameter: --genomes"
    }
    if (!params.gene) {
        error "Missing required parameter: --gene"
    }

    // A value channel lets the single reference genome be reused by many steps.
    ref_ch = Channel
        .fromPath(params.ref, checkIfExists: true)
        .first()

    // Match the original script's sample naming: filename before the first dot.
    genomes_ch = Channel
        .fromPath("${params.genomes}/*", checkIfExists: true)
        .filter { f ->
            def n = f.name.toLowerCase()
            n.endsWith('.fa') || n.endsWith('.fna') || n.endsWith('.fasta')
        }
        .map { genome ->
            def sample = genome.name.split('\\.')[0]
            tuple(sample, genome)
        }

    ANNOTATE_REFERENCE(ref_ch)
    FIND_GENE(ANNOTATE_REFERENCE.out.gff, params.gene)
    EXTRACT_REFERENCE_GENE(ref_ch, FIND_GENE.out.bed, params.gene)

    // Convert the single reference-gene output into a reusable value channel
    // so it can be supplied to every sample-genome task.
    ref_gene_ch = EXTRACT_REFERENCE_GENE.out.fasta.first()

    EXTRACT_SAMPLE_GENE(genomes_ch, ref_gene_ch, params.gene)

    // Annotate the reference gene and every successfully extracted sample gene.
    reference_tuple_ch = ref_gene_ch.map { fasta -> tuple('reference', fasta) }
    all_gene_fastas_ch = EXTRACT_SAMPLE_GENE.out.fasta.mix(reference_tuple_ch)

    ANNOTATE_GENE(all_gene_fastas_ch)

    protein_files_ch = ANNOTATE_GENE.out.protein
        .map { sample, faa -> faa }
        .collect()

    CONCAT_PROTEINS(protein_files_ch)
    ALIGN_PROTEINS(CONCAT_PROTEINS.out.multi)
    CALL_MUTATIONS(ALIGN_PROTEINS.out.aligned, params.gene)
}
