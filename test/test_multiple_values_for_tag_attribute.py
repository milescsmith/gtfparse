import unittest
from io import StringIO

from gtfparse import parse_gtf_and_expand_attributes

class TestMultipleValuesforTagAttribute(unittest.TestCase):
    # failing example from https://github.com/openvax/gtfparse/issues/2
    def setUp(self):
        self.gtf_text = (
            f'1\tprotein_coding\texon\t860260\t860328\t.\t+\t.\t'
            f'gene_id "ENSG00000187634"; transcript_id "ENST00000420190"; '
            f'exon_number "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; '
            f'gene_biotype "protein_coding"; transcript_name "SAMD11-011"; '
            f'transcript_source "havana"; exon_id "ENSE00001637883"; '
            f'tag "cds_end_NF"; tag "mRNA_end_NF"; '
        )


    def test_parse_tag_attributes(self):
        parsed = parse_gtf_and_expand_attributes(StringIO(self.gtf_text))
        tag_column = parsed["tag"]
        self.assertEqual(len(tag_column), 1)
        tags = tag_column[0]
        self.assertEqual(tags, "cds_end_NF,mRNA_end_NF")


    def test_parse_tag_attributes_with_usecols(self):
        parsed = parse_gtf_and_expand_attributes(
            StringIO(self.gtf_text), restrict_attribute_columns=["tag"]
        )
        tag_column = parsed["tag"]
        self.assertEqual(len(tag_column), 1)
        tags = tag_column[0]
        self.assertEqual(tags, "cds_end_NF,mRNA_end_NF")


    def test_parse_tag_attributes_with_usecols_other_column(self):
        parsed = parse_gtf_and_expand_attributes(
            StringIO(self.gtf_text), restrict_attribute_columns=["exon_id"]
        )
        tag_column = parsed.get("tag")

        self.assertIsNone(tag_column, msg=f"Expected 'tag' to get dropped but got {(parsed,)}")

if __name__ == "__main__":
    unittest.main()