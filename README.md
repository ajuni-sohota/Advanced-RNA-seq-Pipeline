# Advanced RNA-Seq Pipeline (ARSP)

A comprehensive, scalable, and modular RNA-Seq analysis pipeline with an interactive Streamlit dashboard for visualization and exploration of results.

## Features

- **Complete RNA-Seq Workflow**: From raw reads to functional analysis
- **Flexible Configuration**: Easily adapt to different experimental designs
- **Multiple Tool Support**: Choose between different aligners, quantifiers, and analysis methods
- **Comprehensive Quality Control**: Monitor quality at every step
- **Advanced Visualization**: Interactive dashboard for exploring results
- **Reproducibility**: Detailed logging and parameter tracking

## Pipeline Workflow


1. **Quality Control**: Raw read quality assessment with FastQC and MultiQC
2. **Read Preprocessing**: Adapter removal and quality trimming with Trimmomatic or fastp
3. **Alignment**: Read alignment with STAR or HISAT2
4. **Quantification**: Gene-level and transcript-level quantification
5. **Expression Analysis**: Normalization, PCA, sample correlation
6. **Differential Expression**: DESeq2 or edgeR analysis
7. **Functional Analysis**: GO term and KEGG pathway enrichment
8. **Reporting**: Comprehensive HTML reports and visualizations

## Interactive Dashboard


The pipeline includes a Streamlit-based interactive dashboard for exploring and visualizing results:

- **Project Overview**: Key statistics and completion status
- **Quality Control**: QC metrics and alignment statistics
- **Expression Analysis**: PCA, sample correlation, and gene expression patterns
- **Differential Expression**: Volcano plots, MA plots, and DE gene tables
- **Functional Enrichment**: GO and KEGG pathway enrichment visualization
- **Export Options**: Download tables and figures for publication

## Installation

### Prerequisites

- Python 3.8 or higher
- R 4.0 or higher
- Required bioinformatics tools:
  - STAR or HISAT2
  - Samtools
  - FastQC
  - MultiQC
  - Trimmomatic or fastp
  - featureCounts (Subread package)
  - Optional: Salmon, Kallisto

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/rnaseq-pipeline.git
   cd rnaseq-pipeline
   ```

2. Create a conda environment (recommended):
   ```bash
   conda env create -f environment.yml
   conda activate rnaseq
   ```

3. Install required R packages:
   ```bash
   Rscript install_r_packages.R
   ```

## Usage

### Running the Pipeline

1. Create a configuration file (see `config_template.yaml` for an example)
2. Run the pipeline:
   ```bash
   python rnaseq_pipeline.py -c config.yaml
   ```

3. To run specific steps:
   ```bash
   python rnaseq_pipeline.py -c config.yaml -s quality_control alignment
   ```

### Running the Dashboard

```bash
streamlit run dashboard.py
```

Then enter the path to your RNA-seq results directory in the dashboard.

## Configuration

The pipeline is controlled by a YAML configuration file with the following sections:

- **project**: General project information
- **samples**: Sample information and metadata
- **references**: Reference genome and annotation files
- **contrasts**: Differential expression comparisons
- **parameters**: Analysis parameters and thresholds
- **tools**: Paths to external tools

See the `config_template.yaml` file for a detailed example.

## Key Considerations

The pipeline incorporates several important considerations for RNA-seq analysis:

- **Batch Effects**: Includes batch in the design matrix for DE analysis
- **Technical Replicates**: Options for handling technical replicates
- **Sequencing Depth Normalization**: Multiple normalization methods (CPM, TPM)
- **Read Strandedness**: Support for unstranded, stranded, and reversely stranded libraries
- **Multi-mapping Reads**: Options for handling multi-mapping reads
- **Count Thresholds**: Filtering low-expressed genes
- **Statistical Rigor**: Multiple testing correction and significance thresholds
- **Experimental Design**: Support for complex designs with multiple factors

## Advanced Features

- **Custom Workflows**: Extend the pipeline with your own analysis modules
- **Parallel Execution**: Multi-threading support for faster analysis
- **Resource Management**: Memory and CPU usage monitoring
- **Error Handling**: Robust error detection and reporting
- **Checkpointing**: Resume from any step after interruption
- **Pipeline Validation**: Built-in checks for input data and parameter consistency

## Output Structure

```
output_dir/
├── 01_fastqc/               # Raw read QC results
├── 02_trimmed/              # Trimmed reads and QC
├── 03_alignment/            # Aligned BAM files and statistics
├── 04_counts/               # Gene and transcript count matrices
├── 05_expression/           # Normalized expression and PCA
├── 06_differential_expression/  # DE analysis results
│   └── functional/         # GO and KEGG enrichment
├── 07_visualization/        # Plots and charts
├── 08_reports/              # Summary reports
└── logs/                    # Detailed pipeline logs
```

## Publication-Ready Figures

The pipeline generates publication-quality figures:

- PCA plots
- Sample correlation heatmaps
- Volcano plots
- MA plots
- Expression heatmaps
- GO enrichment bar plots
- KEGG pathway visualizations

## Citations

If you use this pipeline in your research, please cite the following tools:

- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **STAR**: Dobin A, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21.
- **featureCounts**: Liao Y, et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.
- **DESeq2**: Love MI, et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.
- **clusterProfiler**: Yu G, et al. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5), 284-287.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or issues, please open an issue on GitHub or contact:

Your Name: your.email@example.com

## Acknowledgments

- All the developers of the underlying bioinformatics tools
- The Streamlit team for their amazing visualization library
- The RNA-seq community for guidance on best practices
