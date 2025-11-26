# 🧬 PointMap
A lightweight pipeline for extracting genes and mapping point mutations across bacterial genomes.

PointMap extracts a target gene from multiple genomes, translates the sequences, aligns all proteins, and identifies amino-acid point mutations relative to a user-provided reference gene. It works with any bacterial species and any target gene.

---

## 🚀 Features
- Extract target genes using BLAST + bedtools
- Clean FASTA headers for consistency
- Annotate extracted sequences with Prokka
- Auto-generate the reference protein from the user-provided reference genome or CDS
- Extract amino-acid sequences for all isolates
- Build a combined multi-FASTA protein file
- Align sequences using MAFFT
- Detect amino-acid point mutations relative to the reference
- Fully reproducible using a minimal Conda environment

---

## 📦 Installation

Clone the repository:
```
git clone https://github.com/sbenvari/PointMap.git
cd PointMap
```

Create the environment:
```
conda env create -f environment.yml
conda activate pointmap
```

---

## 📁 Required Inputs

### 1. Directory of genome assemblies  
Supported formats: `.fasta`, `.fa`, `.fna`

Example (included in repo):
```
example_genome/
 └── GCA_000009865.fna
```

### 2. Reference genome or reference gene FASTA  
Must contain the **full CDS** of the target gene.

Example (included in repo):
```
example_ref/
 └── haemo_reference.fna
```

### 3. Gene name  
Used for naming outputs (e.g., `gyrA`, `parC`, `rpoB`).

### 4. Output directory name  
A folder that PointMap will create.

---

## ▶️ Usage

General syntax:
```
./pointmap.sh <GENOMES_DIR> <REFERENCE_FNA> <GENE_NAME> <OUTPUT_DIR>
```

Example — Extract and map mutations in *gyrA*:
```
./pointmap.sh example_genome/ example_ref/haemo_reference.fna gyrA gyrA_results
```

Example — Another gene, same reference:
```
./pointmap.sh example_genome/ example_ref/haemo_reference.fna parC parC_results
```

---

## 📂 Output Structure

```
output_dir/
 ├── 01_sequences/       # Extracted DNA sequences of the target gene
 ├── 02_prokka/          # Prokka annotations of extracted sequences
 ├── 03_proteins/        # Protein sequences (samples + reference)
 ├── 04_alignment/
 │     ├── all.faa       # Combined protein sequences
 │     └── aligned.faa   # MAFFT alignment
 └── 05_results/
       └── <gene>_mutations.txt
```

Example mutation output:
```
isolate1   S84L,E88K
isolate2   No mutations
isolate3   S84A
```

---

## ❗ Reference FASTA Requirements
- Must contain the **full CDS** of the target gene  
- Partial fragments will break mutation numbering  
- A full genome is acceptable (PointMap extracts the gene automatically)  
- The reference is translated and used as the baseline sequence  

---

## 📧 Contact
For issues, bugs, or feature requests, please open an issue on GitHub.
