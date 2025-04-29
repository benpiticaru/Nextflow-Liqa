# Nextflow-Liqa

## Overview
**Nextflow-Liqa** is a bioinformatics pipeline built using [Nextflow](https://www.nextflow.io/) that automates the [LIQA](https://github.com/WGLab/LIQA) long read alternative splicing analysis software. The pipeline includes steps for basecalling, read preprocessing, alignment, isoform quantification, and differential alternative splicing analysis. It is designed to work with high-throughput sequencing data and supports containerized execution using Docker or Singularity.

## Features
- Basecalling using Guppy
- Read preprocessing with Porechop
- Alignment and sorting with Minimap2 and Samtools
- Isoform quantification using LIQA
- Differential expression analysis
- Quality control metrics generation
- GO (Gene Ontology) analysis

## Workflow
The pipeline is implemented as a Nextflow script (`liqa.nf`) and consists of the following steps:
1. **Basecalling**: Converts raw signal data into nucleotide sequences using Guppy.
2. **Read Preprocessing**: Trims adapters and filters low-quality reads using Porechop.
3. **Alignment**: Aligns reads to a reference genome using Minimap2 and sorts the alignments with Samtools.
4. **Isoform Quantification**: Quantifies isoforms and alternative splicing events using LIQA.
5. **Differential Analysis**: Identifies differentially expressed isoforms and splicing events.
6. **Quality Control**: Generates quality control metrics and visualizations.
7. **GO Analysis**: Performs Gene Ontology enrichment analysis on the results.

## Requirements
- **Nextflow**: `>=22.10.0`
- **Docker** or **Singularity** for containerized execution
- **Python**: `>=3.6` with required libraries (e.g., `pandas`, `matplotlib`, `seaborn`)
- **Guppy**, **Minimap2**, **Samtools**, **LIQA**

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/Nextflow-Liqa.git
   cd Nextflow-Liqa
   ```

## Usage
Run the pipeline with the following command:
```bash
nextflow run liqa.nf --input <input_directory> --output <output_directory> --genome <reference_genome>
```

### Parameters
- `--input`: Path to the directory containing raw sequencing data.
- `--output`: Path to the directory where results will be saved.
- `--genome`: Path to the reference genome file.

For additional options, refer to the inline documentation in `liqa.nf`.