# Polars-Bio Agent Skill

![Polars Bio](https://biodatageeks.org/polars-bio/assets/images/logo.png)

## Overview

This repository contains a specialized agent skill for processing large-scale genomic data using **polars-bio** and **Polars**. It demonstrates efficient handling of VCF files, genomic interval joins, and performance benchmarking.

## 📂 Key Scripts

### 1. Analysis & Benchmarking
*   **`scripts/advanced_demo.py`**: 
    *   **Function**: Runs a full analysis pipeline (Load -> Clean -> Join -> Report).
    *   **Features**: Compares **Eager** (in-memory) vs **Streaming** (lazy/out-of-core) execution. Generates `ANALYSIS_REPORT.md` and `clinvar_analysis.html`.
    *   **Usage**: `python3 scripts/advanced_demo.py`

### 2. Utilities (New Implementations)
*   **`scripts/convert_to_parquet.py`**:
    *   **Function**: efficiently converts VCF to Parquet format using streaming. Parquet is much faster for subsequent queries.
    *   **Usage**: `python3 scripts/convert_to_parquet.py`
    
*   **`scripts/fasta_analysis.py`**:
    *   **Function**: Analyzes FASTA files (e.g., Reference Genome) to calculate metrics like GC Content. Demonstrates `scan_fasta`.
    *   **Usage**: `python3 scripts/fasta_analysis.py`

*   **`scripts/annotate_variants.py`**:
    *   **Function**: Annotates VCF variants with external databases (e.g., gnomAD allele frequencies). Uses efficient joining on `chrom, pos, ref, alt`.
    *   **Usage**: `python3 scripts/annotate_variants.py` (Automatically generates a mock gnomAD dataset for demonstration).

### 3. Basics
*   **`scripts/process_vcf.py`**: Simple script for loading and basic interval overlap.

## 📊 Performance Benchmark

The skill includes a rigorous benchmark to demonstrate the benefits of Polars' Streaming engine.

**Results (ClinVar VCF ~177MB):**
*   **Eager Mode**: Loads entire file into RAM. Peak Memory: ~7.8 GB.
*   **Streaming Mode**: Uses `scan_vcf` for lazy loading. Peak Memory: **~5.2 GB** (~33% Reduction).

*See [ANALYSIS_REPORT.md](ANALYSIS_REPORT.md) for detailed charts and metrics.*

## 🚀 Getting Started

1.  **Install Dependencies:**
    ```bash
    pip install polars polars-bio plotly matplotlib seaborn pandas psutil
    ```

2.  **Run the Benchmark:**
    ```bash
    python3 scripts/advanced_demo.py
    ```

3.  **Try Utilities:**
    ```bash
    python3 scripts/convert_to_parquet.py
    python3 scripts/fasta_analysis.py
    ```

## 🔗 Links
*   [Polars-Bio Documentation](https://biodatageeks.org/polars-bio/)
*   [ClinVar Dataset](https://www.ncbi.nlm.nih.gov/clinvar/)