# 🧬 PointMap
A lightweight pipeline for extracting genes and mapping point mutations across bacterial genomes.

PointMap extracts a target gene from multiple genomes, translates the sequences, aligns all proteins, and identifies amino-acid point mutations relative to a user-provided reference gene. It works with any bacterial species and any target gene.

---


## 🚀 Features
- Locate the target gene in both the reference genome and all sample genomes  
  (Prokka annotation + bedtools extraction + BLAST search)
- Extract nucleotide sequences and generate amino-acid sequences for all genomes  
  (Prokka-based translation)
- Align all protein sequences to the reference using MAFFT
- Detect amino-acid point mutations relative to the reference sequence


---

## 📦 Installation

Clone the repository:
```
git clone https://github.com/sbenvari/PointMap.git
cd PointMap
```

Create the environment:
```
conda create -n pointmap
conda activate pointmap
conda env update -f environment.yml
```

---

## 📁 Required Inputs

### 1. Reference genome or reference gene FASTA  
Must contain the **full CDS** of the target gene.

Example (included in repo):
```
example_ref/
 └── haemo_reference.fna
```

### 2. Directory of genome assemblies  
Supported formats: `.fasta`, `.fa`, `.fna`

Example (included in repo):
```
example_genome/
 └── GCA_000009865.fna
```

### 3. Gene name  
Used for naming outputs (e.g., `gyrA`, `parC`, `rpoB`).

### 4. Output directory name  
A folder that PointMap will create.

---

## ▶️ Usage

General syntax:
```
./pointmap.sh <REFERENCE_GENOME> <GENOMES_DIR> <GENE_NAME> <OUTPUT_DIR>
```

Example — Extract and map mutations in *gyrA*:
```
./pointmap.sh path/to/haemo_reference.fna path/to/genomes/ gyrB test_output
```



---

## 📂 Output Structure

```
output_dir/
 ├── sequences_ref/         # Extracted reference gene sequence
 ├── sequences_samples/     # Extracted sample gene sequences
 ├── prokka_out/            # Prokka annotations
 ├── protein/               # Translated protein FASTAs
 ├── multi_sequences.faa    # All proteins concatenated
 ├── aligned_sequences.faa  # MAFFT alignment
 └── mutations_<gene>.txt   # Final mutation list
```

Example mutation output:
```
isolate1   S84L,E88K
isolate2   No mutations
isolate3   S84A
```

---

## 📧 Contact 
For issues, bugs, or feature requests, please open an issue on GitHub.
