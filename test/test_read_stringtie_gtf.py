import unittest
from os.path import exists
from pkg_resources import resource_filename

from gtfparse import read_gtf


class TestRefseqGTF(unittest.TestCase):
    def setUp(self):
        self.b16_gtf_path = resource_filename("test", "data/B16.stringtie.head.gtf")


    def _check_required_columns(self, gtf_df):
        self.assertIn("feature", gtf_df.columns, "Expected column named 'feature' in StringTie GTF")
        self.assertIn("cov", gtf_df.columns, "Expected column named 'cov' in StringTie GTF")
        self.assertIn("FPKM", gtf_df.columns, "Expected column named 'FPKM' in StringTie GTF")
        features = gtf_df["feature"].unique()
        self.assertIn("exon", features, f"No exons in GTF (available: {features})")
        self.assertIn("transcript", features, f"No transcripts in GTF (available: {features})")


    def _check_string_cov_and_FPKM(self, gtf_df):
        for i, feature_name in enumerate(gtf_df["feature"]):
            cov = gtf_df["cov"][i]
            fpkm = gtf_df["FPKM"][i]
            if feature_name == "exon":
                self.assertGreaterEqual(float(cov), 0, msg=f"Expected non-negative cov for exon, got {cov}")
            elif feature_name == "transcript":
                self.assertGreaterEqual(float(cov), 0, msg=f"Expected non-negative cov for transcript, got {cov}")
                self.assertGreaterEqual(float(fpkm), 0, msg=f"Expected non-negative FPKM for transcript, got {fpkm}")


    def _check_float_cov_and_FPKM(self, gtf_df):
        for i, feature_name in enumerate(gtf_df["feature"]):
            cov = gtf_df["cov"][i]
            fpkm = gtf_df["FPKM"][i]
            self.assertIsInstance(cov, float, msg=f"Expected cov to be float but got {cov} : {type(cov)}")
            
            if feature_name == "exon":
                self.assertGreaterEqual(cov, 0, msg=f"Expected non-negative cov for exon, got {cov}")
            elif feature_name == "transcript":
                self.assertIsInstance(fpkm, float, f"Expected FPKM to be float but got {fpkm} : {type(fpkm)}")
                self.assertGreaterEqual(cov, 0, msg=f"Expected non-negative cov for transcript, got {cov}")
                self.assertGreaterEqual(fpkm, 0, msg=f"Expected non-negative FPKM for transcript, got {fpkm}")


    def test_read_stringtie_gtf_as_dataframe(self):
        gtf_df = read_gtf(self.b16_gtf_path, expand_attribute_column=True)
        self._check_required_columns(gtf_df)
        self._check_string_cov_and_FPKM(gtf_df)


    def test_read_stringtie_gtf_as_dataframe_float_values(self):
        gtf_df = read_gtf(self.b16_gtf_path, expand_attribute_column=True, column_converters={"cov": float, "FPKM": float})
        self._check_required_columns(gtf_df)
        self._check_float_cov_and_FPKM(gtf_df)

if __name__ == "__main__":
    unittest.main()