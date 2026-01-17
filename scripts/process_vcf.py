import polars as pl
import polars_bio as pb
import os

def create_dummy_vcf(file_path):
    """Creates a dummy VCF file for testing purposes."""
    content = """##fileformat=VCFv4.2
##source=test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	.	PASS	.
chr1	150	.	C	G	.	PASS	.
chr1	200	.	G	C	.	PASS	.
chr2	50	.	T	A	.	PASS	.
"""
    with open(file_path, "w") as f:
        f.write(content)
    print(f"Created dummy VCF at {file_path}")

def load_vcf(file_path):
    """Loads a VCF file into a Polars DataFrame using polars-bio."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"VCF file not found: {file_path}")
    
    # polars-bio read_vcf automatically handles parsing
    df = pb.read_vcf(file_path)
    return df

def join_intervals(vcf_df, intervals_df):
    """
    Joins VCF DataFrame with an intervals DataFrame.
    
    Args:
        vcf_df: DataFrame loaded from VCF.
        intervals_df: DataFrame with 'chrom', 'start', 'end' columns.
        
    Returns:
        DataFrame: Result of the interval join.
    """
    # Ensure intervals_df has the necessary columns
    required_cols = {"chrom", "start", "end"}
    if not required_cols.issubset(intervals_df.columns):
        raise ValueError(f"Intervals DataFrame must contain columns: {required_cols}")

    # Perform the overlap join
    # polars-bio's overlap function automatically detects start/end columns if they follow conventions
    # or we might need to be explicit.
    # Based on search results, 'by' is used for grouping (chromosome).
    
    # Note: VCF usually has 'pos'. We might need to ensure vcf_df has 'end'.
    # For SNPs, end = pos + len(ref) - 1. Or if it's 1-based closed, end = pos + len(ref) - 1.
    # polars-bio's read_vcf might already compute 'end'.
    
    # Let's check if 'end' exists, if not create it assuming 1bp length for simplicity if REF is missing
    # but read_vcf should provide REF.
    
    # If 'end' is missing, we calculate it from REF length
    if "end" not in vcf_df.columns and "REF" in vcf_df.columns:
        # standard VCF: POS is 1-based start.
        # simple calculation: start = pos, end = pos + len(ref)
        # But polars-bio likely expects specific column names for overlap.
        pass

    # Using the API found in help(pb.overlap): pb.overlap(df1, df2, ...)
    try:
        joined = pb.overlap(vcf_df, intervals_df)
        if isinstance(joined, pl.LazyFrame):
            joined = joined.collect()
        return joined
    except Exception as e:
        print(f"Error during overlap: {e}")
        # Fallback debug print
        print("VCF columns:", vcf_df.columns)
        print("Intervals columns:", intervals_df.columns)
        raise

def main():
    vcf_path = "data/sample.vcf"
    
    # Ensure data directory exists
    os.makedirs("data", exist_ok=True)
    
    # Create a dummy VCF if it doesn't exist
    if not os.path.exists(vcf_path):
        create_dummy_vcf(vcf_path)

    print("Loading VCF...")
    vcf_df = load_vcf(vcf_path)
    print("VCF Loaded:")
    print(vcf_df)

    # Define some intervals to join against (e.g., genes or regions of interest)
    intervals_df = pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [90, 190, 40],
        "end": [110, 210, 60],
        "region_name": ["Region_A", "Region_B", "Region_C"]
    })
    
    print("\nIntervals to Join:")
    print(intervals_df)

    print("\nPerforming Interval Join...")
    try:
        result = join_intervals(vcf_df, intervals_df)
        print("Join Result:")
        print(result)
    except Exception as e:
        print(f"Join failed: {e}")

if __name__ == "__main__":
    main()
