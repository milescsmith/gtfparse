import unittest
from io import StringIO

import numpy as np
import pandas as pd

from gtfparse import (REQUIRED_COLUMNS, ParsingError, parse_gtf,
                      parse_gtf_and_expand_attributes)


class TestParseGTF(unittest.TestCase):
    def setUp(self):
        self.gtf_text = (
            "# sample GTF data copied from:\n"
            "# http://useast.ensembl.org/info/website/upload/gff.html?redirect=no\n"
            "1\ttranscribed_unprocessed_pseudogene\tgene\t11869\t"
            '14409\t.\t+\t.\tgene_id "ENSG00000223972"; '
            'gene_name "DDX11L1"; gene_source "havana"; gene_biotype '
            '"transcribed_unprocessed_pseudogene";\n'
            "1\tprocessed_transcript\ttranscript\t11869\t"
            '14409\t.\t+\t.\tgene_id "ENSG00000223972";'
            'transcript_id "ENST00000456328"; gene_name "DDX11L1";'
            'gene_source "havana";'
            'gene_biotype "transcribed_unprocessed_pseudogene";'
            'transcript_name "DDX11L1-002";transcript_source "havana";'
        )
        self.expanded_columns = pd.Index(
            REQUIRED_COLUMNS[:8]
            + [
                "gene_id",
                "gene_name",
                "gene_source",
                "gene_biotype",
                "transcript_id",
                "transcript_name",
                "transcript_source",
            ]
        )
        self.required_columns = pd.Index(REQUIRED_COLUMNS)
        self.expected_seqname = pd.Series(["1", "1"], dtype=str)
        self.expected_start = pd.Series([11869, 11869], dtype=np.int32)
        self.expected_end = pd.Series([14409, 14409], dtype=np.int32)
        self.expected_gene_id = pd.Series(
            ["ENSG00000223972", "ENSG00000223972"], dtype=object
        )
        self.expected_transcript_id = pd.Series(["", "ENST00000456328"], dtype=object)
        self.expected_attribute = pd.Series(
            [
                'gene_id "ENSG00000223972";gene_name "DDX11L1";gene_source "havana";gene_biotype "transcribed_unprocessed_pseudogene";',
                'gene_id "ENSG00000223972";transcript_id "ENST00000456328";gene_name "DDX11L1";gene_source "havana";gene_biotype "transcribed_unprocessed_pseudogene";transcript_name "DDX11L1-002";transcript_source "havana";',
            ],
            dtype=str,
        )

    def test_parse_gtf_lines_with_expanded_attributes(self):
        parsed_df = parse_gtf_and_expand_attributes(StringIO(self.gtf_text))

        self.assertTrue(
            parsed_df.columns.equals(self.expanded_columns),
            msg=f"columns not equal. Expected {self.expanded_columns} got {parsed_df.columns}",
        )
        self.assertTrue(
            parsed_df["seqname"].equals(self.expected_seqname),
            msg=f"transcript_id not equal. Expected {self.expected_seqname} got {parsed_df['seqname']}",
        )

        self.assertTrue(
            parsed_df["start"].equals(self.expected_start),
            msg=f"transcript_id not equal. Expected {self.expected_start} got {parsed_df['start']}",
        )
        self.assertTrue(
            parsed_df["end"].equals(self.expected_end),
            msg=f"transcript_id not equal. Expected {self.expected_end} got {parsed_df['end']}",
        )

        self.assertTrue(np.isnan(parsed_df["score"]).all(), msg="Unexpected scores")
        self.assertTrue(
            parsed_df["gene_id"].equals(self.expected_gene_id),
            msg=f"gene_id not equal. Expected {self.expected_gene_id} got {parsed_df['gene_id']}",
        )
        self.assertTrue(
            parsed_df["transcript_id"].equals(self.expected_transcript_id),
            msg=f"transcript_id not equal. Expected {self.expected_transcript_id} got {parsed_df['transcript_id']}",
        )

    def test_parse_gtf_lines_without_expand_attributes(self):
        parsed_df = parse_gtf(StringIO(self.gtf_text))

        self.assertTrue(parsed_df.columns.equals(self.required_columns))
        self.assertTrue(parsed_df["seqname"].equals(self.expected_seqname))

        self.assertTrue(parsed_df["start"].equals(self.expected_start))
        self.assertTrue(parsed_df["end"].equals(self.expected_end))

        self.assertTrue(np.isnan(parsed_df["score"]).all())
        self.assertTrue(parsed_df["attribute"].equals(self.expected_attribute))

    def test_parse_gtf_lines_error_too_many_fields(self):
        bad_gtf_text = self.gtf_text.replace(" ", "\t")
        # pylint: disable=no-value-for-parameter
        with self.assertRaises(ParsingError):
            parse_gtf(StringIO(bad_gtf_text))

    def test_parse_gtf_lines_error_too_few_fields(self):
        bad_gtf_text = self.gtf_text.replace("\t", " ")
        # pylint: disable=no-value-for-parameter
        with self.assertRaises(ParsingError):
            parse_gtf(StringIO(bad_gtf_text))


if __name__ == "__main__":
    unittest.main()
