from io import StringIO

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from gtfparse.parsing_error import ParsingError
from gtfparse.read_gtf import parse_gtf, parse_gtf_and_expand_attributes
from gtfparse.required_columns import REQUIRED_COLUMNS

# ruff: noqa: S101


@pytest.fixture
def gtf_text() -> str:
    return (
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


@pytest.fixture
def parsed_gtf(gtf_text: str) -> pd.DataFrame:
    return parse_gtf(StringIO(gtf_text))


@pytest.fixture
def expanded_parsed_gtf(gtf_text: str) -> pd.DataFrame:
    return parse_gtf_and_expand_attributes(StringIO(gtf_text))


@pytest.fixture
def expanded_columns() -> pd.Index:
    return pd.Index(
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


@pytest.fixture
def required_columns() -> pd.Index:
    return pd.Index(REQUIRED_COLUMNS)


@pytest.fixture
def expected_seqname() -> pd.Series:
    return pd.Series(["1", "1"], name="seqname", dtype=str)


@pytest.fixture
def expected_start() -> pd.Series:
    return pd.Series([11869, 11869], name="start", dtype=np.int64)


@pytest.fixture
def expected_end() -> pd.Series:
    return pd.Series([14409, 14409], name="end", dtype=np.int64)


@pytest.fixture
def expected_score() -> pd.Series:
    return pd.Series([".", "."], name="score", dtype=object)


@pytest.fixture
def expected_gene_id() -> pd.Series:
    return pd.Series(["ENSG00000223972", "ENSG00000223972"], name="gene_id", dtype=object)


@pytest.fixture
def expected_transcript_id() -> pd.Series:
    return pd.Series(["", "ENST00000456328"], name="transcript_id", dtype=object)


@pytest.fixture
def expected_attribute() -> pd.Series:
    return pd.Series(
        [
            'gene_id "ENSG00000223972";gene_name "DDX11L1";gene_source "havana";gene_biotype "transcribed_unprocessed_pseudogene";', # noqa: E501
            'gene_id "ENSG00000223972";transcript_id "ENST00000456328";gene_name "DDX11L1";gene_source "havana";gene_biotype "transcribed_unprocessed_pseudogene";transcript_name "DDX11L1-002";transcript_source "havana";', # noqa: E501
        ],
        name="attribute",
        dtype=str,
    )


def test_expanded_parsed_gtf_columns(expanded_parsed_gtf: pd.DataFrame, expanded_columns: pd.Index):
    pdt.assert_index_equal(expanded_parsed_gtf.columns, expanded_columns)


def test_expanded_parsed_gtf_seqname(expanded_parsed_gtf: pd.DataFrame, expected_seqname: pd.Series):
    pdt.assert_series_equal(expanded_parsed_gtf["seqname"], expected_seqname)


def test_expanded_parsed_gtf_start(expanded_parsed_gtf: pd.DataFrame, expected_start: pd.Series):
    pdt.assert_series_equal(expanded_parsed_gtf["start"], expected_start)


def test_expanded_parsed_gtf_end(expanded_parsed_gtf: pd.DataFrame, expected_end: pd.Series):
    pdt.assert_series_equal(expanded_parsed_gtf["end"], expected_end)


def test_expanded_parsed_gtf_score(expanded_parsed_gtf: pd.DataFrame, expected_score: pd.Series):
    pdt.assert_series_equal(expanded_parsed_gtf["score"], expected_score)


def test_expanded_parsed_gtf_gene_id(expanded_parsed_gtf: pd.DataFrame, expected_gene_id: pd.Series):
    pdt.assert_series_equal(expanded_parsed_gtf["gene_id"], expected_gene_id)


def test_expanded_parsed_gtf_transcript_id(expanded_parsed_gtf: pd.DataFrame, expected_transcript_id: pd.Series):
    pdt.assert_series_equal(expanded_parsed_gtf["transcript_id"], expected_transcript_id)


def test_parsed_gtf_columns(parsed_gtf: pd.DataFrame, required_columns: pd.Index):
    pdt.assert_index_equal(parsed_gtf.columns, required_columns)


def test_parsed_gtf_seqname(parsed_gtf: pd.DataFrame, expected_seqname: pd.Series):
    pdt.assert_series_equal(parsed_gtf["seqname"], expected_seqname)


def test_parsed_gtf_start(parsed_gtf: pd.DataFrame, expected_start: pd.Series):
    pdt.assert_series_equal(parsed_gtf["start"], expected_start)


def test_parsed_gtf_end(parsed_gtf: pd.DataFrame, expected_end: pd.Series):
    pdt.assert_series_equal(parsed_gtf["end"], expected_end)


def test_parsed_gtf_score(parsed_gtf: pd.DataFrame, expected_score: pd.Series):
    pdt.assert_series_equal(parsed_gtf["score"], expected_score)


def test_parsed_gtf_attribute(parsed_gtf: pd.DataFrame, expected_attribute: pd.Series):
    pdt.assert_series_equal(parsed_gtf["attribute"], expected_attribute)


@pytest.fixture
def bad_gtf_text_too_many_fields(gtf_text: str) -> str:
    return gtf_text.replace(" ", "\t")


def test_parse_bad_gtf_error_too_many_fields(bad_gtf_text_too_many_fields: str):
    with pytest.raises(ParsingError):
        parse_gtf(StringIO(bad_gtf_text_too_many_fields))


@pytest.fixture
def bad_gtf_text_too_few_fields(gtf_text: str) -> str:
    return gtf_text.replace("\t", " ")


def test_parse_bad_gtf_error_too_few_fields(bad_gtf_text_too_few_fields: str):
    with pytest.raises(ParsingError):
        parse_gtf(StringIO(bad_gtf_text_too_few_fields))
