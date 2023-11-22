from io import StringIO

import pandas as pd
import pytest

from gtfparse.read_gtf import parse_gtf_and_expand_attributes


# ruff: noqa: S101
# failing example from https://github.com/openvax/gtfparse/issues/2
@pytest.fixture
def gtf_text() -> str:
    return (
        "1\tprotein_coding\texon\t860260\t860328\t.\t+\t.\t"
        'gene_id "ENSG00000187634"; transcript_id "ENST00000420190"; '
        'exon_number "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; '
        'gene_biotype "protein_coding"; transcript_name "SAMD11-011"; '
        'transcript_source "havana"; exon_id "ENSE00001637883"; '
        'tag "cds_end_NF"; tag "mRNA_end_NF"; '
    )


@pytest.fixture
def parsed_gtf(gtf_text: str) -> pd.DataFrame:
    return parse_gtf_and_expand_attributes(StringIO(gtf_text))


def test_parse_tag_attributes(parsed_gtf: pd.DataFrame):
    assert len(parsed_gtf["tag"]) == 1
    assert parsed_gtf["tag"][0] == "cds_end_NF,mRNA_end_NF"


@pytest.fixture
def parsed_restricted_tag_col_gtf(gtf_text: str) -> pd.DataFrame:
    return parse_gtf_and_expand_attributes(StringIO(gtf_text), restrict_attribute_columns=["tag"])


def test_parse_tag_attributes_with_usecols(parsed_restricted_tag_col_gtf: pd.DataFrame):
    assert len(parsed_restricted_tag_col_gtf["tag"]) == 1
    assert parsed_restricted_tag_col_gtf["tag"][0] == "cds_end_NF,mRNA_end_NF"


@pytest.fixture
def parsed_restricted_exon_id_col_gtf(gtf_text: str) -> pd.DataFrame:
    return parse_gtf_and_expand_attributes(StringIO(gtf_text), restrict_attribute_columns=["exon_id"])


def test_parse_tag_attributes_with_usecols_other_column(parsed_restricted_exon_id_col_gtf: pd.DataFrame):
    assert parsed_restricted_exon_id_col_gtf.get("tag") is None
