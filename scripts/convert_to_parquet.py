import polars as pl
import polars_bio as pb
import time
import os

def convert_vcf_to_parquet(vcf_path, parquet_path):
    """
    Converts a VCF file to Parquet format using streaming.
    """
    print(f"Converting {vcf_path} to {parquet_path}...")
    t0 = time.time()
    
    # scan_vcf is lazy and memory efficient
    lf = pb.scan_vcf(vcf_path)
    
    # sink_parquet processes in chunks/streams
    lf.sink_parquet(parquet_path)
    
    dt = time.time() - t0
    print(f"Conversion complete in {dt:.2f}s")
    
    # Verify result
    df = pl.read_parquet(parquet_path)
    print(f"Parquet Shape: {df.height} rows, {df.width} columns")
    print(f"Parquet Size: {os.path.getsize(parquet_path) / (1024*1024):.2f} MB")

def main():
    # Use existing sample data
    vcf_path = "polars-bio-agent-skill/data/clinvar.vcf.gz"
    parquet_path = "polars-bio-agent-skill/data/clinvar.parquet"
    
    if not os.path.exists(vcf_path):
        print(f"Error: {vcf_path} not found. Please run advanced_demo.py first or download data.")
        return

    convert_vcf_to_parquet(vcf_path, parquet_path)

if __name__ == "__main__":
    main()
