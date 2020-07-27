import unittest
from os.path import exists
from pkg_resources import resource_filename

from gtfparse import read_gtf

REFSEQ_GTF_PATH = data_path("refseq.ucsc.small.gtf")


class TestRefseGTF(unittest.TestCase):
    def setUp(self):
        self.refseq = resource_filename("data", "refseq.ucsc.small.gtf")

        self.assertTrue(exists(self.refseq), msg="Cannot find 'refseq.ucsc.small.gtf'")

    def _check_required_columns(self):
        refseq = read_gtf(self.refseq)

        self.assertIn("feature", refseq.columns, msg="Expected column named 'feature' not found in RefSeq GTF")
        self.assertIn("gene_id", refseq.columns, msg="Expected column named 'gene_id' not found in RefSeq GTF")
        self.assertIn("transcript_id", refseq.columns, msg="Expected column named 'transcript_id' not found in RefSeq GTF")
        
        self.assertIn("exon", refseq['feature'].unique(), msg=f"No exon features in GTF (available: {refseq['feature'].unique()})")
        self.assertIn("CDS", refseq['feature'].unique(), msg=f"No CDS features in GTF (available: {refseq['feature'].unique()})")
