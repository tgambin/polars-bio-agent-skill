import polars as pl
import polars_bio as pb
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os

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
        # Ensure images directory exists
        os.makedirs("polars-bio-agent-skill/images", exist_ok=True)
        path = f"polars-bio-agent-skill/images/{filename}"
        plt.savefig(path, bbox_inches='tight', dpi=100)
        plt.close()
        return f"images/{filename}"
    except Exception as e:
        print(f"Failed to save plot {filename}: {e}")
        return None

def main():
    start_total = time.time()
    print("Starting Advanced Polars-Bio Demo...")
    
    # Paths
    vcf_path = "polars-bio-agent-skill/data/clinvar.vcf.gz"
    bed_path = "polars-bio-agent-skill/data/cytoBand.txt.gz"
    output_html = "polars-bio-agent-skill/clinvar_analysis.html"

    # --- 1. Data Loading ---
    print("\n[1/4] Loading Data...")
    t0 = time.time()
    try:
        vcf_df = pb.read_vcf(vcf_path)
    except Exception as e:
        print(f"Failed to load VCF: {e}")
        return

    load_time_vcf = time.time() - t0
    print(f"  - VCF Loaded in {load_time_vcf:.2f}s")
    print(f"  - VCF Shape: {vcf_df.height} rows, {vcf_df.width} columns")
    
    # Load Cytobands
    t0 = time.time()
    cytoband_df = pl.read_csv(
        bed_path, 
        has_header=False, 
        separator="\t",
        new_columns=["chrom", "start", "end", "band", "stain"]
    )
    load_time_bed = time.time() - t0
    print(f"  - Cytobands Loaded in {load_time_bed:.2f}s")
    print(f"  - Cytobands Shape: {cytoband_df.height} rows")

    # --- 2. Preprocessing & Feature Engineering ---
    print("\n[2/4] Preprocessing & Feature Engineering...")
    t0 = time.time()
    
    cols = vcf_df.columns
    clnsig_col = next((c for c in cols if c.upper() == "CLNSIG"), None)
    info_col = next((c for c in cols if c.upper() == "INFO"), None)

    if clnsig_col:
        vcf_df = vcf_df.with_columns(pl.col(clnsig_col).alias("CLNSIG_Simple"))
    elif info_col:
        vcf_df = vcf_df.with_columns(
            pl.col(info_col).str.extract(r"CLNSIG=([^;]+)", 1).alias("CLNSIG_Simple")
        )
    else:
        vcf_df = vcf_df.with_columns(pl.lit("Unknown").alias("CLNSIG_Simple"))

    if vcf_df.schema["CLNSIG_Simple"] == pl.List(pl.Utf8):
        vcf_df = vcf_df.with_columns(
            pl.col("CLNSIG_Simple").list.get(0).alias("CLNSIG_Simple")
        )

    vcf_df = vcf_df.with_columns(
        pl.col("CLNSIG_Simple").str.split("/").list.get(0).str.split("|").list.get(0).alias("clnsig_simple")
    )
    
    pathogenic_df = vcf_df.filter(
        pl.col("clnsig_simple").str.to_lowercase().str.contains("pathogenic") &
        ~pl.col("clnsig_simple").str.to_lowercase().str.contains("conflicting")
    )
    
    pathogenic_df = pathogenic_df.select(["chrom", "start", "end", "clnsig_simple"])
    cytoband_df = cytoband_df.select(["chrom", "start", "end", "band"])
    
    prep_time = time.time() - t0
    print(f"  - Preprocessing finished in {prep_time:.2f}s")
    print(f"  - Pathogenic Variants Found: {pathogenic_df.height}")

    if pathogenic_df.height == 0:
        print("  - No pathogenic variants found. Skipping Join and Visualization.")
        return

    # --- 3. Interval Join (Overlap) ---
    print("\n[3/4] Performing Interval Join (Variants x Cytobands)...")
    t0 = time.time()
    
    vcf_chr = str(vcf_df["chrom"][0])
    cyto_chr = str(cytoband_df["chrom"][0])
    
    if not vcf_chr.startswith("chr") and cyto_chr.startswith("chr"):
        vcf_df = vcf_df.with_columns(("chr" + pl.col("chrom")).alias("chrom"))
        pathogenic_df = pathogenic_df.with_columns(("chr" + pl.col("chrom")).alias("chrom"))
    elif vcf_chr.startswith("chr") and not cyto_chr.startswith("chr"):
        cytoband_df = cytoband_df.with_columns(("chr" + pl.col("chrom")).alias("chrom"))

    joined_df = pb.overlap(pathogenic_df, cytoband_df)
    
    if isinstance(joined_df, pl.LazyFrame):
        joined_df = joined_df.collect()

    join_time = time.time() - t0
    print(f"  - Interval Join finished in {join_time:.2f}s")
    print(f"  - Joined Rows: {joined_df.height}")

    # --- 4. Visualization & Reporting ---
    print("\n[4/4] Generating Visualizations...")
    t0 = time.time()
    
    chrom_counts = vcf_df.group_by("chrom").len().sort("chrom")
    clnsig_counts = vcf_df.group_by("clnsig_simple").len().sort("len", descending=True).head(15)
    
    cols = joined_df.columns
    band_col = next((c for c in cols if "band" in c), None)
    chrom_col = "chrom_1" if "chrom_1" in cols else "chrom"
    
    band_counts = pl.DataFrame()
    if band_col:
        print(f"  - Using '{band_col}' as band column.")
        band_counts = joined_df.group_by([chrom_col, band_col]).len().sort("len", descending=True).head(20)
        band_counts = band_counts.with_columns(
            (pl.col(chrom_col) + pl.col(band_col)).alias("FullBand")
        )

    # 4a. Interactive HTML (Plotly)
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "xy"}, {"type": "domain"}],
               [{"type": "xy", "colspan": 2}, None]],
        subplot_titles=("Variants per Chromosome", "Clinical Significance (Top 15)", "Top 20 Cytobands with Pathogenic Variants")
    )
    chrom_counts_p = chrom_counts.to_pandas()
    fig.add_trace(go.Bar(x=chrom_counts_p["chrom"], y=chrom_counts_p["len"], name="Variants"), row=1, col=1)
    
    clnsig_p = clnsig_counts.to_pandas()
    fig.add_trace(go.Pie(labels=clnsig_p["clnsig_simple"], values=clnsig_p["len"], name="Clinical Sig"), row=1, col=2)
    
    band_p = None
    if not band_counts.is_empty():
        band_p = band_counts.to_pandas()
        fig.add_trace(go.Bar(x=band_p["FullBand"], y=band_p["len"], name="Pathogenic Vars"), row=2, col=1)
        fig.update_xaxes(title_text="Cytoband", row=2, col=1)

    fig.update_layout(title_text="ClinVar VCF Analysis with Polars-Bio", height=800, showlegend=True)
    fig.write_html(output_html)
    
    # 4b. Static Plots (Matplotlib/Seaborn)
    print("  - Generating static PNG plots for Markdown...")
    sns.set_theme(style="whitegrid")
    
    def plot_chroms():
        plt.figure(figsize=(12, 6))
        sns.barplot(data=chrom_counts_p, x="chrom", y="len", color="steelblue")
        plt.title("Variants per Chromosome")
        plt.xticks(rotation=45)
        plt.tight_layout()
    img_chrom = save_static_plot(plot_chroms, "chrom_counts.png")

    def plot_clnsig():
        plt.figure(figsize=(10, 6))
        sns.barplot(data=clnsig_p, y="clnsig_simple", x="len", orient="h", palette="viridis")
        plt.title("Clinical Significance Distribution (Top 15)")
        plt.xlabel("Count")
        plt.tight_layout()
    img_clnsig = save_static_plot(plot_clnsig, "clnsig_dist.png")

    img_bands = None
    if not band_counts.is_empty():
        def plot_bands():
            plt.figure(figsize=(12, 8))
            sns.barplot(data=band_p, x="len", y="FullBand", orient="h", color="coral")
            plt.title("Top 20 Cytobands with Pathogenic Variants")
            plt.xlabel("Count")
            plt.tight_layout()
        img_bands = save_static_plot(plot_bands, "top_bands.png")

    viz_time = time.time() - t0
    print(f"  - Visualization generated in {viz_time:.2f}s")
    
    # --- 5. Generate Markdown Report ---
    md_path = "polars-bio-agent-skill/ANALYSIS_REPORT.md"
    print(f"\n[5/5] Generating Markdown Report ({md_path})...")
    
    with open(md_path, "w") as f:
        f.write("# ClinVar VCF Analysis Report\n\n")
        f.write(f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 1. Dataset Overview\n")
        f.write(f"- **Total Variants in VCF:** {vcf_df.height:,}\n")
        f.write(f"- **Pathogenic Variants:** {pathogenic_df.height:,}\n")
        f.write(f"- **Joined with Cytobands:** {joined_df.height:,}\n\n")
        
        if img_chrom:
            f.write(f"![Variants per Chromosome]({img_chrom})\n\n")
        
        f.write("## 2. Clinical Significance Distribution (Top 15)\n")
        if img_clnsig:
             f.write(f"![Clinical Significance]({img_clnsig})\n\n")
             
        f.write("| Clinical Significance | Count |\n")
        f.write("| :--- | :--- |\n")
        for row in clnsig_counts.iter_rows():
            f.write(f"| {row[0]} | {row[1]:,} |\n")
        f.write("\n")
        
        if not band_counts.is_empty():
            f.write("## 3. Top 20 Cytobands with Pathogenic Variants\n")
            f.write(f"*> Based on interval overlap with UCSC Cytobands*\n\n")
            
            if img_bands:
                f.write(f"![Top Cytobands]({img_bands})\n\n")
                
            f.write("| Rank | Chromosome | Band | Full Name | Pathogenic Variant Count |\n")
            f.write("| :--- | :--- | :--- | :--- | :--- |\n")
            
            for i, row in enumerate(band_counts.iter_rows(named=True), 1):
                chrom = row[chrom_col]
                band = row[band_col]
                count = row["len"]
                full = row["FullBand"]
                f.write(f"| {i} | {chrom} | {band} | {full} | {count:,} |\n")
        
        f.write("\n---\n")
        f.write("*Note: For interactive visualizations, download `clinvar_analysis.html`.*")

    print(f"  - Markdown report saved to: {md_path}")

    total_time = time.time() - start_total
    print(f"\nTotal Execution Time: {total_time:.2f}s")
    
    if os.path.exists(output_html):
        print(f"HTML Output size: {format_bytes(os.path.getsize(output_html))}")

if __name__ == "__main__":
    main()