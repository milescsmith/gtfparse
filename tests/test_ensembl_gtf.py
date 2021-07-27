import unittest

# import pandas as pd
import modin.pandas as pd
from pkg_resources import resource_filename

from gtfparse import read_gtf


class TestEnsemblGTF(unittest.TestCase):
    def setUp(self):
        self.ensembl_gtf = resource_filename("tests", "data/ensembl_grch37.head.gtf")
        self.expected_features = pd.Series(
            ["gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon"],
            dtype=str,
        )
        # first 1000 lines of GTF only contained these genes
        self.expected_gene_names = pd.Series(
            [
                "FAM41C",
                "CICP27",
                "RNU6-1100P",
                "NOC2L",
                "AP006222.1",
                "LINC01128",
                "RP4-669L17.1",
                "RP11-206L10.2",
                "PLEKHN1",
                "WBP1LP7",
                "RP5-857K21.1",
                "RP5-857K21.5",
                "RNU6-1199P",
                "RP11-206L10.10",
                "RP11-54O7.16",
                "CICP7",
                "AL627309.1",
                "RP5-857K21.11",
                "DDX11L1",
                "RP5-857K21.3",
                "RP11-34P13.7",
                "AL669831.1",
                "MTATP6P1",
                "CICP3",
                "WBP1LP6",
                "LINC00115",
                "hsa-mir-6723",
                "RP5-857K21.7",
                "SAMD11",
                "RP11-206L10.5",
                "RP11-34P13.8",
                "RP11-206L10.9",
                "RP11-34P13.15",
                "TUBB8P11",
                "MTATP8P1",
                "RP4-669L17.8",
                "RP11-206L10.1",
                "RP11-34P13.13",
                "RP11-206L10.3",
                "RP11-206L10.4",
                "RP11-54O7.3",
                "RP5-857K21.2",
                "OR4F5",
                "MTND1P23",
                "AL645608.1",
                "RP11-34P13.16",
                "RP11-34P13.14",
                "AP006222.2",
                "OR4F29",
                "RP4-669L17.4",
                "AL732372.1",
                "OR4G4P",
                "MTND2P28",
                "OR4F16",
                "KLHL17",
                "FAM138A",
                "OR4G11P",
                "FAM87B",
                "RP5-857K21.15",
                "AL645608.2",
                "RP11-206L10.8",
                "RP5-857K21.4",
                "MIR1302-10",
                "RP11-54O7.2",
                "RP4-669L17.10",
                "RP11-54O7.1",
                "RP11-34P13.9",
                "WASH7P",
                "RP4-669L17.2",
            ],
            dtype=str,
        )

    def test_ensembl_gtf_columns(self):
        df = read_gtf(self.ensembl_gtf, expand_attribute_column=True)
        features = pd.Series(df["feature"].str.strip('"').unique().astype(str))
        self.assertTrue(
            all(features.isin(self.expected_features)),
            msg="Unexpected features found in ensembl_gtf",
        )

    def test_ensembl_gtf_gene_names(self):
        df = read_gtf(self.ensembl_gtf, expand_attribute_column=True)
        gene_names = pd.Series(df["gene_name"].str.strip('"').unique().astype(str))
        self.assertTrue(
            all(gene_names.isin(self.expected_gene_names)),
            msg=(
                f"Wrong gene names: {gene_names}, "
                f"missing {self.expected_gene_names[~self.expected_gene_names.isin(gene_names)]} "
                f"and unexpected {gene_names[~gene_names.isin(self.expected_gene_names)]}"
            ),
        )

    def test_ensembl_gtf_gene_names_with_usecols(self):
        df = read_gtf(
            self.ensembl_gtf, usecols=["gene_name"], expand_attribute_column=True
        )
        gene_names = pd.Series(df["gene_name"].str.strip('"').unique().astype(str))
        self.assertTrue(
            all(gene_names.isin(self.expected_gene_names)),
            msg=(
                f"Wrong gene names: {gene_names}, "
                f"missing {self.expected_gene_names[~self.expected_gene_names.isin(gene_names)]} "
                f"and unexpected {gene_names[~gene_names.isin(self.expected_gene_names)]}"
            ),
        )

    def test_ensembl_gtf_gene_names_zip(self):
        df = read_gtf(self.ensembl_gtf + ".gz", expand_attribute_column=True)
        gene_names = pd.Series(df["gene_name"].str.strip('"').unique().astype(str))
        self.assertTrue(
            all(gene_names.isin(self.expected_gene_names)),
            msg=(
                f"Wrong gene names: {gene_names}, "
                f"missing {self.expected_gene_names[~self.expected_gene_names.isin(gene_names)]} "
                f"and unexpected {gene_names[~gene_names.isin(self.expected_gene_names)]}"
            ),
        )

    def test_ensembl_gtf_gene_names_with_usecols_gzip(self):
        df = read_gtf(
            self.ensembl_gtf + ".gz",
            usecols=["gene_name"],
            expand_attribute_column=True,
        )
        gene_names = pd.Series(df["gene_name"].str.strip('"').unique().astype(str))
        self.assertTrue(
            all(gene_names.isin(self.expected_gene_names)),
            msg=(
                f"Wrong gene names: {gene_names}, "
                f"missing {self.expected_gene_names[~self.expected_gene_names.isin(gene_names)]} "
                f"and unexpected {gene_names[~gene_names.isin(self.expected_gene_names)]}"
            ),
        )
