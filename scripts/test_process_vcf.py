import unittest
import polars as pl
import polars_bio as pb
import os
import shutil
import tempfile
from process_vcf import load_vcf, join_intervals, create_dummy_vcf

class TestProcessVcf(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.dummy_vcf_path = os.path.join(self.test_dir, "test.vcf")
        create_dummy_vcf(self.dummy_vcf_path)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_load_vcf(self):
        df = load_vcf(self.dummy_vcf_path)
        self.assertIsInstance(df, pl.DataFrame)
        self.assertEqual(df.height, 4)
        self.assertIn("chrom", df.columns)
        self.assertIn("start", df.columns)
        # Check content of first row
        self.assertEqual(df["start"][0], 100)

    def test_join_intervals(self):
        vcf_df = load_vcf(self.dummy_vcf_path)
        
        intervals_df = pl.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [90, 190],
            "end": [110, 210],
            "region_name": ["Region_A", "Region_B"]
        })
        
        result = join_intervals(vcf_df, intervals_df)
        
        if isinstance(result, pl.LazyFrame):
            result = result.collect()
            
        self.assertEqual(result.height, 2)
        
        # Check that we found the expected overlaps
        starts = []
        if "start_1" in result.columns:
            starts.extend(result["start_1"].to_list())
        if "start_2" in result.columns:
            starts.extend(result["start_2"].to_list())
            
        self.assertIn(100, starts)
        self.assertIn(200, starts)
        self.assertNotIn(150, starts)

    def test_join_intervals_no_overlap(self):
        vcf_df = load_vcf(self.dummy_vcf_path)
        
        intervals_df = pl.DataFrame({
            "chrom": ["chr1"],
            "start": [5000],
            "end": [6000],
            "region_name": ["Far_Away"]
        })
        
        result = join_intervals(vcf_df, intervals_df)
        if isinstance(result, pl.LazyFrame):
            result = result.collect()
            
        self.assertEqual(result.height, 0)

if __name__ == "__main__":
    unittest.main()
