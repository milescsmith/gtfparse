from importlib.resources import as_file, files

import pandas as pd
import pytest

from gtfparse.read_gtf import read_gtf

# ruff: noqa: S101

@pytest.fixture
def refseq():
    with as_file(files("tests.data").joinpath("refseq.ucsc.small.gtf")) as gtf:
        return read_gtf(gtf, expand_attribute_column=True)


def test_refseq_columns(refseq):
    assert (
        pd.Series(
            [
                "feature",
                "gene_id",
                "transcript_id",
            ]
        )
        .isin(refseq.columns)
        .all()
    )


def test_refseq_features(refseq):
    assert (
        pd.Series(
            [
                "exon",
                "CDS",
            ]
        )
        .isin(refseq["feature"])
        .all
    )
