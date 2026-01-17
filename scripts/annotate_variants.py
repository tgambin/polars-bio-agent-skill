import polars as pl
import polars_bio as pb
import numpy as np
import os
import time

def create_mock_gnomad(vcf_path, out_path):
    print("Creating mock gnomAD dataset from ClinVar sample...")
    if not os.path.exists(vcf_path):
         raise FileNotFoundError(f"VCF not found: {vcf_path}")
         
    vcf = pb.read_vcf(vcf_path).select(["chrom", "start", "ref", "alt"]).unique()
    mock_db = vcf.sample(n=50000, seed=42)
    af_values = np.random.beta(a=0.5, b=5, size=mock_db.height)
    mock_db = mock_db.with_columns(pl.Series("gnomad_af", af_values))
    mock_db.write_parquet(out_path)

def annotate_variants(vcf_path, db_path):
    print("\nAnnotating ClinVar with gnomAD (Mock)...")
    t0 = time.time()
    
    vcf = pb.scan_vcf(vcf_path)
    gnomad = pl.scan_parquet(db_path)
    
    annotated = vcf.join(
        gnomad,
        on=["chrom", "start", "ref", "alt"],
        how="left"
    )
    
    # Calculate stats
    stats = annotated.select([
        pl.len().alias("total"),
        pl.col("gnomad_af").is_not_null().sum().alias("annotated_count")
    ]).collect()
    
    dt = time.time() - t0
    
    total = stats["total"][0]
    annotated_count = stats["annotated_count"][0]
    
    return {
        "Time (s)": dt,
        "Total Variants": total,
        "Annotated Variants": annotated_count,
        "Annotation Rate": annotated_count / total if total > 0 else 0
    }

def main():
    vcf_path = "polars-bio-agent-skill/data/clinvar.vcf.gz"
    db_path = "polars-bio-agent-skill/data/gnomad_mock.parquet"
    
    if not os.path.exists(vcf_path):
        print("VCF not found.")
        return

    # Create mock DB if not exists
    if not os.path.exists(db_path):
        create_mock_gnomad(vcf_path, db_path)
        
    metrics = annotate_variants(vcf_path, db_path)
    print(metrics)

if __name__ == "__main__":
    main()
