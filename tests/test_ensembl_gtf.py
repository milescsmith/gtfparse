from importlib.resources import as_file, files

import pandas as pd
import pytest
from loguru import logger

from gtfparse.logging import init_logger
from gtfparse.read_gtf import read_gtf

# ruff: noqa: S101
logger.disable("gtfparse")
init_logger(verbose=3, save=False)
logger.enable("gtfparse")


@pytest.fixture
def ensembl_gtf() -> pd.DataFrame:
    with as_file(files("tests.data").joinpath("ensembl_grch37.head.gtf")) as gtf:
        return read_gtf(gtf, expand_attribute_column=True)


@pytest.fixture
def gzipped_ensembl_gtf() -> pd.DataFrame:
    with as_file(files("tests.data").joinpath("ensembl_grch37.head.gtf.gz")) as gtf:
        return read_gtf(gtf, expand_attribute_column=True)


@pytest.fixture
def gzipped_usecol_ensembl_gtf() -> pd.DataFrame:
    with as_file(files("tests.data").joinpath("ensembl_grch37.head.gtf.gz")) as gtf:
        return read_gtf(gtf, usecols=["gene_name"], expand_attribute_column=True)


@pytest.fixture
def expected_features() -> pd.Series:
    return pd.Series(
        ["gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon"],
        dtype=str,
    )


@pytest.fixture
def expected_gene_names() -> pd.Series:
    return pd.Series(
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


# ensembl_gtf == df
@pytest.fixture
def features(ensembl_gtf: pd.DataFrame):
    return pd.Series(ensembl_gtf["feature"].str.strip('"').unique().astype(str))


def test_ensembl_gtf_columns(features: pd.Series, expected_features: pd.Series):
    assert all(features.isin(expected_features))


def test_ensembl_gtf_gene_names(ensembl_gtf: pd.DataFrame, expected_gene_names: pd.Series):
    gene_names = pd.Series(ensembl_gtf["gene_name"].str.strip('"').unique().astype(str))
    assert all(gene_names.isin(expected_gene_names))


def test_ensembl_gtf_gene_names_with_usecols(ensembl_gtf: pd.DataFrame, expected_gene_names: pd.Series):
    gene_names = pd.Series(ensembl_gtf["gene_name"].str.strip('"').unique().astype(str))
    assert all(gene_names.isin(expected_gene_names))


def test_ensembl_gtf_gene_names_zip(gzipped_ensembl_gtf: pd.DataFrame, expected_gene_names: pd.Series):
    gene_names = pd.Series(gzipped_ensembl_gtf["gene_name"].str.strip('"').unique().astype(str))
    assert all(gene_names.isin(expected_gene_names))


def test_ensembl_gtf_gene_names_with_usecols_gzip(
    gzipped_usecol_ensembl_gtf: pd.DataFrame, expected_gene_names: pd.Series
):
    gene_names = pd.Series(gzipped_usecol_ensembl_gtf["gene_name"].str.strip('"').unique().astype(str))
    assert all(gene_names.isin(expected_gene_names))
