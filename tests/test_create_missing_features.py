import unittest
from io import StringIO

from gtfparse import create_missing_features, parse_gtf_and_expand_attributes

# two lines from the Ensembl 54 human GTF containing only a stop_codon and
# exon features, but from which gene and transcript information could be
# inferred


class TestCreateMissingFeatures(unittest.TestCase):
    def setUp(self):
        self.gtf_text = (
            f'# seqname biotype feature start end score strand frame '
            f'attribute\n'
            f'18\tprotein_coding\tstop_codon\t32630766\t'
            f'32630768\t.\t-\t0\t'
            f'gene_id "ENSG00000134779"; transcript_id "ENST00000334295"; '
            f'exon_number "7";gene_name "C18orf10";transcript_name '
            f'"C18orf10-201"\n'
            f'18\tprotein_coding\texon\t32663078\t32663157\t'
            f'.\t+\t.\tgene_id "ENSG00000150477"; '
            f'transcript_id "ENST00000383055"; exon_number "1"; gene_name '
            f'"KIAA1328"; transcript_name "KIAA1328-202";'
        )

        self.gtf_df = parse_gtf_and_expand_attributes(StringIO(self.gtf_text))

    def test_create_missing_features_identity(self):
        df_should_be_same = create_missing_features(self.gtf_df, {})
        self.assertEqual(
            len(self.gtf_df),
            len(df_should_be_same),
            msg="GTF DataFrames are not the same size",
        )

    def _check_expanded_dataframe(self, df):
        self.assertIn(
            "gene",
            df["feature"].unique(),
            msg="Extended GTF should contain gene feature",
        )
        self.assertIn(
            "transcript",
            df["feature"].unique(),
            msg="Extended GTF should contain transcript feature",
        )

        C18orf10_201_transcript_mask = (df["feature"] == "transcript") & (
            df["transcript_name"] == 'C18orf10-201'
        )
        self.assertEqual(
            len(df[C18orf10_201_transcript_mask]),
            1,
            msg=f"Expected only 1 gene entry for C18orf10-201, got {df[C18orf10_201_transcript_mask]}",
        )

        transcript_seqname = df[C18orf10_201_transcript_mask].seqname.iloc[0]
        self.assertEqual(
            transcript_seqname,
            "18",
            msg=f"Wrong seqname for C18orf10-201: {transcript_seqname}",
        )

        transcript_start = df[C18orf10_201_transcript_mask].start.iloc[0]
        self.assertEqual(
            transcript_start,
            32630766,
            msg=f"Wrong start for C18orf10-201: {transcript_start}",
        )

        transcript_end = df[C18orf10_201_transcript_mask].end.iloc[0]
        self.assertEqual(
            transcript_end,
            32630768,
            msg=f"Wrong end for C18orf10-201: {transcript_end}",
        )

        transcript_strand = df[C18orf10_201_transcript_mask].strand.iloc[0]
        self.assertEqual(
            transcript_strand,
            "-",
            msg=f"Wrong strand for C18orf10-201: {transcript_strand}",
        )

        KIAA1328_gene_mask = (df["feature"] == "gene") & (df["gene_name"] == 'KIAA1328')
        self.assertEqual(
            len(df[KIAA1328_gene_mask]),
            1,
            msg=f"Expected only 1 gene entry for KIAA1328, found {len(df[KIAA1328_gene_mask])}",
        )

        gene_seqname = df[KIAA1328_gene_mask].seqname.iloc[0]
        self.assertEqual(
            gene_seqname, "18", msg=f"Wrong seqname for KIAA1328: {gene_seqname}"
        )

        gene_start = df[KIAA1328_gene_mask].start.iloc[0]
        self.assertEqual(
            gene_start, 32663078, msg=f"Wrong start for KIAA1328: {(gene_start,)}"
        )

        gene_end = df[KIAA1328_gene_mask].end.iloc[0]
        self.assertEqual(
            gene_end, 32663157, msg=f"Wrong end for KIAA1328: {(gene_end,)}"
        )

        gene_strand = df[KIAA1328_gene_mask].strand.iloc[0]
        self.assertEqual(
            gene_strand, "+", msg=f"Wrong strand for KIAA1328: {gene_strand}"
        )

    def test_create_missing_features(self):
        self.assertNotIn(
            "gene",
            self.gtf_df["feature"].unique(),
            msg="Original GTF should not contain gene feature",
        )
        self.assertNotIn(
            "transcript",
            self.gtf_df["feature"].unique(),
            msg="Original GTF should not contain transcript feature",
        )

        df_extra_features = create_missing_features(
            self.gtf_df,
            unique_keys={"gene": "gene_id", "transcript": "transcript_id"},
            extra_columns={
                "gene": {"gene_name"},
                "transcript": {"gene_id", "gene_name", "transcript_name"},
            },
        )
        self._check_expanded_dataframe(df_extra_features)


if __name__ == "__main__":
    unittest.main()
