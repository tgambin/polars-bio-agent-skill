# Polars-Bio Agent Skill

![Polars Bio](https://biodatageeks.org/polars-bio/assets/images/logo.png)

## Overview

This repository contains a specialized agent skill for processing large-scale genomic data using **polars-bio** and **Polars**. It demonstrates efficient handling of VCF files (Variant Call Format) and performing genomic interval joins (overlaps) with annotations like cytobands.

## 🎯 Goal

The primary goal of this skill is to provide a robust, high-performance toolkit for:
*   Loading and parsing VCF files (~100MB+).
*   Cleaning and feature engineering genomic data.
*   Performing interval overlap joins with other genomic regions (BED files).
*   Benchmarking execution modes (Eager vs. Streaming) to optimize performance.

## 📂 Key Files

*   **[SKILL.md](SKILL.md)**: Detailed definition of the skill's capabilities and quick start guide.
*   **[ANALYSIS_REPORT.md](ANALYSIS_REPORT.md)**: A generated report containing analysis results, visualizations, and performance benchmarks.
*   **`scripts/advanced_demo.py`**: The main Python script that runs the pipeline, performs the benchmark, and generates the report.
*   **`scripts/process_vcf.py`**: A simpler, modular script for basic VCF processing.

## 🚀 Usage

### Prerequisites

Ensure you have the required dependencies installed:

```bash
pip install polars polars-bio plotly matplotlib seaborn pandas
```

### Running the Benchmark & Demo

To run the full analysis pipeline, including the performance benchmark (Eager vs. Streaming) and report generation:

```bash
python3 scripts/advanced_demo.py
```

This will:
1.  Load the ClinVar VCF (~177MB) and UCSC Cytobands.
2.  Run the processing pipeline in **Eager** mode.
3.  Run the processing pipeline in **Streaming** mode.
4.  Generate interactive charts in `clinvar_analysis.html`.
5.  Generate a static Markdown report with plots in `ANALYSIS_REPORT.md`.

### Running Tests

To verify the core functionality:

```bash
python3 scripts/test_process_vcf.py
```

## 📊 Performance Benchmark

The `advanced_demo.py` script compares two execution strategies:

1.  **Eager Mode**: Materializes DataFrames at each step (Load -> Filter -> Join). Good for debugging and smaller data.
2.  **Streaming Mode**: Uses Polars' lazy evaluation and streaming engine (`collect(streaming=True)`). Designed for datasets larger than RAM.

**View the full results and charts in the [Analysis Report](ANALYSIS_REPORT.md).**

## 🔗 Links

*   [Polars-Bio Documentation](https://biodatageeks.org/polars-bio/)
*   [Polars Documentation](https://pola.rs/)
*   [ClinVar Dataset](https://www.ncbi.nlm.nih.gov/clinvar/)
