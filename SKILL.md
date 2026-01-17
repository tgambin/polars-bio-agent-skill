---
name: polars-bio
description: A skill for processing genomic data using polars-bio, focusing on VCF loading and interval joins.
license: Apache-2.0
---

# Polars-Bio Skill

## Overview

This skill demonstrates how to use `polars-bio` (and `polars`) to efficiently process genomic data. It specifically covers loading Variant Call Format (VCF) files and performing interval joins (overlaps) with other genomic intervals.

## Capabilities

- **VCF Ingestion:** High-performance loading of VCF files into Polars DataFrames.
- **Interval Joins:** Efficient overlap joins between genomic intervals (e.g., finding which variants fall within specific genomic regions).

## Quick Start

```python
import polars as pl
import polars_bio as pb

# Load a VCF file
vcf_df = pb.read_vcf("path/to/file.vcf")

# Create or load intervals (e.g., from a BED file or DataFrame)
regions_df = pl.DataFrame({
    "chrom": ["chr1", "chr1"],
    "start": [100, 200],
    "end": [150, 250],
    "region_name": ["region_1", "region_2"]
})

# Perform an interval join (overlap)
# Finds variants in vcf_df that overlap with regions_df
joined_df = vcf_df.pb.overlap(
    regions_df,
    on=["chrom", "start", "end"], # Columns to join on
    method="inner" # or "left", "outer"
)

print(joined_df)
```

## Dependencies

- `polars`
- `polars-bio`
