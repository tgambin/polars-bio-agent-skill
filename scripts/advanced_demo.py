import polars as pl
import polars_bio as pb
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os
import pandas as pd
import psutil
import threading
import multiprocessing
import sys

# Append scripts directory to path to allow importing utilities
sys.path.append(os.path.join(os.path.dirname(__file__)))
from convert_to_parquet import convert_vcf_to_parquet
from fasta_analysis import analyze_fasta
from annotate_variants import annotate_variants, create_mock_gnomad

def format_bytes(size):
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} TB"

def save_static_plot(fig_func, filename):
    """Helper to save matplotlib figures."""
    try:
        fig_func()
        os.makedirs("polars-bio-agent-skill/images", exist_ok=True)
        path = f"polars-bio-agent-skill/images/{filename}"
        plt.savefig(path, bbox_inches='tight', dpi=100)
        plt.close()
        return f"images/{filename}"
    except Exception as e:
        print(f"Failed to save plot {filename}: {e}")
        return None

def monitor_memory(pid, stop_event, result_list):
    """Monitors RSS memory of a process."""
    process = psutil.Process(pid)
    peak = 0
    while not stop_event.is_set():
        try:
            mem = process.memory_info().rss
            if mem > peak:
                peak = mem
            time.sleep(0.05)
        except psutil.NoSuchProcess:
            break
    result_list.append(peak)

def run_pipeline_task(mode, vcf_path, bed_path, result_queue):
    """
    Task to be run in a separate process.
    Performs the pipeline and puts metrics into the queue.
    """
    # Start Memory Monitor
    mem_peak_container = []
    stop_mem_event = threading.Event()
    mem_thread = threading.Thread(target=monitor_memory, args=(os.getpid(), stop_mem_event, mem_peak_container))
    mem_thread.start()

    try:
        metrics = {"Mode": mode, "Load_Time": 0, "Prep_Time": 0, "Join_Time": 0, "Total_Time": 0}
        t_start_pipeline = time.time()
        
        # --- 1. Load ---
        t0 = time.time()
        
        if mode == 'Streaming':
            # Use scan_vcf for lazy loading
            vcf = pb.scan_vcf(vcf_path)
            # Cytoband is small, but let's be consistent
            cytoband = pl.scan_csv(
                bed_path, has_header=False, separator="\t",
                new_columns=["chrom", "start", "end", "band", "stain"]
            )
        else:
            # Eager Mode
            vcf = pb.read_vcf(vcf_path)
            cytoband = pl.read_csv(
                bed_path, has_header=False, separator="\t",
                new_columns=["chrom", "start", "end", "band", "stain"]
            )
            
            # Ensure Materialization for Eager Benchmarking
            if isinstance(vcf, pl.LazyFrame):
                vcf = vcf.collect()

        metrics["Load_Time"] = time.time() - t0
        
        vcf_shape_str = f"{vcf.height} rows" if isinstance(vcf, pl.DataFrame) else "Lazy"

        # --- 2. Preprocess ---
        t0 = time.time()
        
        cols = vcf.columns
        clnsig_col = next((c for c in cols if c.upper() == "CLNSIG"), None)
        info_col = next((c for c in cols if c.upper() == "INFO"), None)

        if clnsig_col:
            vcf = vcf.with_columns(pl.col(clnsig_col).alias("CLNSIG_Simple"))
        elif info_col:
            vcf = vcf.with_columns(
                pl.col(info_col).str.extract(r"CLNSIG=([^;]+)", 1).alias("CLNSIG_Simple")
            )
        else:
            vcf = vcf.with_columns(pl.lit("Unknown").alias("CLNSIG_Simple"))

        if vcf.schema["CLNSIG_Simple"] == pl.List(pl.Utf8):
            vcf = vcf.with_columns(pl.col("CLNSIG_Simple").list.get(0).alias("CLNSIG_Simple"))

        vcf = vcf.with_columns(
            pl.col("CLNSIG_Simple").str.split("/").list.get(0).str.split("|").list.get(0).alias("clnsig_simple")
        )
        
        pathogenic = vcf.filter(
            pl.col("clnsig_simple").str.to_lowercase().str.contains("pathogenic") &
            ~pl.col("clnsig_simple").str.to_lowercase().str.contains("conflicting")
        )
        
        pathogenic = pathogenic.select(["chrom", "start", "end", "clnsig_simple"])
        cytoband = cytoband.select(["chrom", "start", "end", "band"])
        
        metrics["Prep_Time"] = time.time() - t0

        # --- 3. Join ---
        t0 = time.time()
        
        # Simplified logic for benchmark: assume transform needed
        pathogenic = pathogenic.with_columns(("chr" + pl.col("chrom")).alias("chrom"))
        
        joined = pb.overlap(pathogenic, cytoband)
        
        if mode == 'Eager':
            if isinstance(joined, pl.LazyFrame):
                joined = joined.collect()
        else:
            # Streaming
            # engine='streaming' is now standard, but explicit arg is safe for benchmark distinction
            # if supported by version. If not, use standard collect which optimizes.
            try:
                joined = joined.collect(streaming=True)
            except:
                joined = joined.collect()

        metrics["Join_Time"] = time.time() - t0
        metrics["Total_Time"] = time.time() - t_start_pipeline
        metrics["Result_Rows"] = joined.height
        
        # Stop Memory Monitor
        stop_mem_event.set()
        mem_thread.join()
        
        metrics["Memory_Peak_MB"] = mem_peak_container[0] / (1024 * 1024)
        
        result_queue.put(metrics)
        
    except Exception as e:
        result_queue.put({"error": str(e)})

def run_pipeline_isolated(mode, vcf_path, bed_path):
    """Wraps pipeline execution in a separate process."""
    q = multiprocessing.Queue()
    p = multiprocessing.Process(target=run_pipeline_task, args=(mode, vcf_path, bed_path, q))
    p.start()
    p.join()
    return q.get()

def main():
    print("Starting Advanced Polars-Bio Demo & Benchmark...")
    
    # Paths (relative to root)
    base_dir = "polars-bio-agent-skill"
    vcf_path = os.path.join(base_dir, "data/clinvar.vcf.gz")
    bed_path = os.path.join(base_dir, "data/cytoBand.txt.gz")
    output_html = os.path.join(base_dir, "clinvar_analysis.html")
    
    # Utility paths
    parquet_path = os.path.join(base_dir, "data/clinvar.parquet")
    fasta_path = os.path.join(base_dir, "data/chr22.fa.gz")
    mock_db_path = os.path.join(base_dir, "data/gnomad_mock.parquet")

    results = []
    
    # Run Eager
    print("\n--- Benchmarking Eager Mode ---")
    res_eager = run_pipeline_isolated('Eager', vcf_path, bed_path)
    if "error" in res_eager:
        print(f"Eager failed: {res_eager['error']}")
        return
    results.append(res_eager)
    
    # Run Streaming
    print("\n--- Benchmarking Streaming Mode ---")
    res_stream = run_pipeline_isolated('Streaming', vcf_path, bed_path)
    if "error" in res_stream:
        print(f"Streaming failed: {res_stream['error']}")
        return
    results.append(res_stream)
    
    # Convert results to DataFrame
    bench_df = pd.DataFrame(results)
    
    # --- Run Utilities ---
    print("\n--- Running Utility Demos ---")
    
    # 1. Parquet Conversion
    print("Running Parquet Conversion...")
    parquet_metrics = convert_vcf_to_parquet(vcf_path, parquet_path)
    
    # 2. FASTA Analysis
    print("Running FASTA Analysis...")
    fasta_metrics = analyze_fasta(fasta_path)
    
    # 3. Annotation
    print("Running Variant Annotation...")
    if not os.path.exists(mock_db_path):
        create_mock_gnomad(vcf_path, mock_db_path)
    annotation_metrics = annotate_variants(vcf_path, mock_db_path)

    # --- Visualizations for Report ---
    
    # 1. Execution Time Comparison
    def plot_time():
        plt.figure(figsize=(10, 6))
        melted = bench_df.melt(id_vars="Mode", value_vars=["Load_Time", "Prep_Time", "Join_Time"], var_name="Stage", value_name="Time(s)")
        sns.barplot(data=melted, x="Stage", y="Time(s)", hue="Mode", palette="muted")
        plt.title("Execution Time by Stage")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
    img_time = save_static_plot(plot_time, "benchmark_time.png")

    # 2. Total Time & Memory Footprint
    def plot_memory():
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        # Bar for Time
        sns.barplot(data=bench_df, x="Mode", y="Total_Time", color="skyblue", ax=ax1, alpha=0.6)
        ax1.set_ylabel("Total Time (s)", color="blue")
        ax1.tick_params(axis='y', labelcolor="blue")
        
        # Line/Point for Memory
        ax2 = ax1.twinx()
        sns.lineplot(data=bench_df, x="Mode", y="Memory_Peak_MB", color="green", marker="s", ax=ax2, linewidth=2, markersize=10)
        ax2.set_ylabel("Peak Memory (MB)", color="green")
        ax2.tick_params(axis='y', labelcolor="green")
        ax2.set_ylim(0, max(bench_df["Memory_Peak_MB"]) * 1.2)
        
        plt.title("Performance: Time vs Memory Footprint")
    img_mem = save_static_plot(plot_memory, "benchmark_memory.png")

    # --- Re-run Aggregations for Plots ---
    print("\nReloading data for visualizations...")
    vcf = pb.scan_vcf(vcf_path) # Use scan (streaming logic) for viz data too, faster
    cytoband = pl.scan_csv(bed_path, has_header=False, separator="\t", new_columns=["chrom", "start", "end", "band", "stain"])
    
    cols = vcf.collect_schema().names()
    clnsig_col = next((c for c in cols if c.upper() == "CLNSIG"), None)
    if clnsig_col: vcf = vcf.with_columns(pl.col(clnsig_col).alias("clnsig_simple"))
    
    # Note: Lazy handling of List
    # We apply transform without checking first row if we trust schema
    # But schema check for list needs collect_schema
    schema = vcf.collect_schema()
    if schema.get("clnsig_simple") == pl.List(pl.Utf8):
        vcf = vcf.with_columns(pl.col("clnsig_simple").list.get(0))
        
    vcf = vcf.with_columns(pl.col("clnsig_simple").str.split("/").list.get(0).str.split("|").list.get(0))
    pathogenic = vcf.filter(
         pl.col("clnsig_simple").str.to_lowercase().str.contains("pathogenic") &
        ~pl.col("clnsig_simple").str.to_lowercase().str.contains("conflicting")
    )
    pathogenic = pathogenic.select(["chrom", "start", "end", "clnsig_simple"])
    cytoband = cytoband.select(["chrom", "start", "end", "band"])
    pathogenic = pathogenic.with_columns(("chr" + pl.col("chrom")).alias("chrom"))
    
    df_viz = pb.overlap(pathogenic, cytoband).collect()
    
    # Visualizations
    chrom_col = "chrom_1" if "chrom_1" in df_viz.columns else "chrom"
    band_col = next((c for c in df_viz.columns if "band" in c), "band")
    
    # Plot 1: Chromosomes
    chrom_counts_p = df_viz.group_by(chrom_col).len().sort(chrom_col).to_pandas()
    def plot_chroms():
        plt.figure(figsize=(12, 6))
        sns.barplot(data=chrom_counts_p, x=chrom_col, y="len", color="steelblue")
        plt.title("Pathogenic Variants per Chromosome (Joined)")
        plt.xticks(rotation=45)
    img_chrom = save_static_plot(plot_chroms, "chrom_counts.png")

    # Plot 2: Top Bands
    band_counts = df_viz.group_by([chrom_col, band_col]).len().sort("len", descending=True).head(20)
    band_counts = band_counts.with_columns((pl.col(chrom_col) + pl.col(band_col)).alias("FullBand"))
    band_p = band_counts.to_pandas()
    def plot_bands():
        plt.figure(figsize=(12, 8))
        sns.barplot(data=band_p, x="len", y="FullBand", orient="h", color="coral")
        plt.title("Top 20 Cytobands with Pathogenic Variants")
    img_bands = save_static_plot(plot_bands, "top_bands.png")

    
    # --- Generate Markdown Report ---
    md_path = os.path.join(base_dir, "ANALYSIS_REPORT.md")
    print(f"\nGenerating Markdown Report ({md_path})...")
    
    with open(md_path, "w") as f:
        f.write("# Polars-Bio Performance & Capabilities Report\n\n")
        f.write(f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 1. Executive Summary\n")
        f.write(f"This report evaluates the **polars-bio-agent-skill** on the ClinVar dataset.\n")
        f.write("It includes a performance benchmark of execution modes and demonstrates key utility features.\n\n")
        
        # --- Section 2: Benchmark ---
        f.write("## 2. Performance Benchmark (Eager vs Streaming)\n")
        f.write("### Metrics Comparison\n")
        
        cols = bench_df.columns
        header = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        rows = []
        for _, row in bench_df.iterrows():
            row_vals = []
            for val in row:
                if isinstance(val, float):
                    row_vals.append(f"{val:.2f}")
                else:
                    row_vals.append(str(val))
            rows.append("| " + " | ".join(row_vals) + " |")
        f.write("\n".join([header, sep] + rows))
        f.write("\n\n")
        
        f.write("### Memory Footprint Analysis\n")
        if img_mem: f.write(f"![Memory Footprint]({img_mem})\n\n")
        
        # --- Section 3: Utility Demos ---
        f.write("## 3. Utility Capabilities Demonstration\n")
        f.write("The skill includes specialized tools for common genomic tasks. Below are the execution results:\n\n")
        
        f.write("### A. High-Performance Parquet Conversion\n")
        f.write(f"Converted {vcf_path} to optimal Parquet format.\n\n")
        f.write("| Metric | Value |\n| :--- | :--- |\n")
        for k, v in parquet_metrics.items():
            val = f"{v:.2f}" if isinstance(v, float) else str(v)
            f.write(f"| {k} | {val} |\n")
        f.write("\n")
        
        f.write("### B. FASTA Sequence Analysis\n")
        f.write(f"Analyzed {fasta_path} (GC Content calculation).\n\n")
        f.write("| Metric | Value |\n| :--- | :--- |\n")
        for k, v in fasta_metrics.items():
            val = f"{v:.2f}" if isinstance(v, float) else str(v)
            f.write(f"| {k} | {val} |\n")
        f.write("\n")
        
        f.write("### C. Large-Scale Variant Annotation\n")
        f.write(f"Annotated ClinVar with mock gnomAD allele frequencies.\n\n")
        f.write("| Metric | Value |\n| :--- | :--- |\n")
        for k, v in annotation_metrics.items():
            val = f"{v:.2%}" if "Rate" in k else (f"{v:.2f}" if isinstance(v, float) else f"{v:,}")
            f.write(f"| {k} | {val} |\n")
        f.write("\n")

        # --- Section 4: Genomic Analysis ---
        f.write("## 4. Genomic Analysis (Pathogenic Variants)\n")
        f.write(f"- **Total Overlapping Variants:** {df_viz.height:,}\n\n")
        
        f.write("### Distribution by Chromosome\n")
        if img_chrom: f.write(f"![Chromosomes]({img_chrom})\n\n")
        
        f.write("### Top Cytobands\n")
        if img_bands: f.write(f"![Bands]({img_bands})\n\n")

    print(f"Report saved to: {md_path}")

if __name__ == "__main__":
    main()
