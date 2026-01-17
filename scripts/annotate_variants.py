import polars as pl
import polars_bio as pb
import numpy as np
import os

def create_mock_gnomad(vcf_path, out_path):
    print("Creating mock gnomAD dataset from ClinVar sample...")
    # Read ClinVar to get real loci
    vcf = pb.read_vcf(vcf_path).select(["chrom", "start", "ref", "alt"]).unique()
    
    # Sample 50,000 variants to simulate coverage
    mock_db = vcf.sample(n=50000, seed=42)
    
    # Add random Allele Frequency (AF)
    # Using numpy to generate random float array
    af_values = np.random.beta(a=0.5, b=5, size=mock_db.height)
    
    mock_db = mock_db.with_columns(
        pl.Series("gnomad_af", af_values)
    )
    
    # Save as Parquet
    mock_db.write_parquet(out_path)
    print(f"Mock gnomAD saved to {out_path} ({mock_db.height} variants)")

def annotate_variants(vcf_path, db_path):
    print("\nAnnotating ClinVar with gnomAD (Mock)...")
    
    # Lazy read both
    vcf = pb.scan_vcf(vcf_path)
    gnomad = pl.scan_parquet(db_path)
    
    # VCF uses 'start' for position.
    # Join Keys: chrom, start, ref, alt
    # Note: Ensure types match (chrom usually str)
    
    # Standard Polars Join
    annotated = vcf.join(
        gnomad,
        on=["chrom", "start", "ref", "alt"],
        how="left"
    )
    
    # Filter for annotated variants (where gnomad_af is not null)
    # Collect a sample
    result = annotated.filter(pl.col("gnomad_af").is_not_null()) \
                      .select(["chrom", "start", "ref", "alt", "gnomad_af"])
                      .head(10)
                      .collect()
    
    print("Sample Annotated Variants:")
    print(result)
    
    # Stats
    total = vcf.select(pl.len()).collect().item()
    annotated_count = annotated.filter(pl.col("gnomad_af").is_not_null()).select(pl.len()).collect().item()
    print(f"\nTotal Variants: {total}")
    print(f"Annotated Variants: {annotated_count} ({annotated_count/total:.1%})")

def main():
    vcf_path = "polars-bio-agent-skill/data/clinvar.vcf.gz"
    db_path = "polars-bio-agent-skill/data/gnomad_mock.parquet"
    
    if not os.path.exists(vcf_path):
        print("VCF not found.")
        return

    # Create mock DB if not exists
    if not os.path.exists(db_path):
        create_mock_gnomad(vcf_path, db_path)
        
    annotate_variants(vcf_path, db_path)

if __name__ == "__main__":
    main()
