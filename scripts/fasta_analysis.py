import polars as pl
import polars_bio as pb
import time
import os

def analyze_fasta(fasta_path):
    print(f"Analyzing FASTA: {fasta_path}")
    t0 = time.time()
    
    # Lazy read
    lf = pb.scan_fasta(fasta_path)
    
    # Calculate GC Content
    lf = lf.with_columns(
        pl.col("sequence").str.to_uppercase().alias("seq_upper")
    )
    
    lf = lf.with_columns(
        (pl.col("seq_upper").str.count_matches("G") + pl.col("seq_upper").str.count_matches("C")).alias("gc_count"),
        pl.col("seq_upper").str.len_chars().alias("seq_len")
    )
    
    lf = lf.with_columns(
        (pl.col("gc_count") / pl.col("seq_len")).alias("gc_content")
    )
    
    # Aggregate results (Weighted average GC)
    result = lf.select([
        pl.col("seq_len").sum().alias("total_len"),
        pl.col("gc_count").sum().alias("total_gc")
    ]).collect()
    
    avg_gc = result["total_gc"][0] / result["total_len"][0]
    dt = time.time() - t0
    
    return {
        "Time (s)": dt,
        "Total Sequence Length (bp)": result["total_len"][0],
        "Average GC Content": avg_gc
    }

def main():
    fasta_path = "polars-bio-agent-skill/data/chr22.fa.gz"
    
    if not os.path.exists(fasta_path):
        print(f"Error: {fasta_path} not found.")
        return

    metrics = analyze_fasta(fasta_path)
    print(metrics)

if __name__ == "__main__":
    main()
