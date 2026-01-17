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
    # Schema usually: id, description, seq
    # We calculate (count(G) + count(C)) / len(seq)
    
    # Note: Sequences in FASTA can be UPPER or lower case.
    lf = lf.with_columns(
        pl.col("seq").str.to_uppercase().alias("seq_upper")
    )
    
    lf = lf.with_columns(
        (pl.col("seq_upper").str.count_matches("G") + pl.col("seq_upper").str.count_matches("C")).alias("gc_count"),
        pl.col("seq_upper").str.len_chars().alias("seq_len")
    )
    
    lf = lf.with_columns(
        (pl.col("gc_count") / pl.col("seq_len")).alias("gc_content")
    )
    
    # Aggregation: Average GC content across all sequences (e.g. contigs/chromosomes)
    # For chr22.fa.gz, there might be one main sequence and random contigs.
    
    result = lf.select(["id", "seq_len", "gc_content"]).collect()
    
    dt = time.time() - t0
    print(f"Analysis complete in {dt:.2f}s")
    print(result)

def main():
    fasta_path = "polars-bio-agent-skill/data/chr22.fa.gz"
    
    if not os.path.exists(fasta_path):
        print(f"Error: {fasta_path} not found.")
        return

    analyze_fasta(fasta_path)

if __name__ == "__main__":
    main()
