from importlib.resources import as_file, files
from io import StringIO
from sys import intern
from typing import TextIO

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from loguru import logger

from gtfparse.create_missing_features import create_missing_features
from gtfparse.logging import init_logger
from gtfparse.read_gtf import parse_frame, parse_gtf_and_expand_attributes

# two lines from the Ensembl 54 human GTF containing only a stop_codon and
# exon features, but from which gene and transcript information could be
# inferred

# ruff: noqa: S101

logger.disable("gtfparse")
init_logger(verbose=3, save=False)
logger.enable("gtfparse")


@pytest.fixture
def saved_gtf_df() -> pd.DataFrame:
    gtf_as_df = files("tests.data").joinpath("gtf_as_df.csv")
    with as_file(gtf_as_df) as data:
        return pd.read_csv(
            data,
            index_col=0,
            converters={
                "frame": parse_frame,
                "seqname": intern,
                "source": intern,
                "feature": intern,
                "score": intern,
            },
            dtype={"start": np.int64, "end": np.int64, "exon_number": np.int64, "strand": "category"},
        )


@pytest.fixture
def gtf_text() -> TextIO:
    return StringIO(
        "# seqname biotype feature start end score strand frame "
        "attribute\n"
        "18\tprotein_coding\tstop_codon\t32630766\t"
        "32630768\t.\t-\t0\t"
        'gene_id "ENSG00000134779"; transcript_id "ENST00000334295"; '
        'exon_number "7";gene_name "C18orf10";transcript_name '
        '"C18orf10-201"\n'
        "18\tprotein_coding\texon\t32663078\t32663157\t"
        '.\t+\t.\tgene_id "ENSG00000150477"; '
        'transcript_id "ENST00000383055"; exon_number "1"; gene_name '
        '"KIAA1328"; transcript_name "KIAA1328-202";'
    )


@pytest.fixture
def created_gtf_df(gtf_text) -> pd.DataFrame:
    return parse_gtf_and_expand_attributes(gtf_text)


def test_parse_gtf_and_expand_attributes(saved_gtf_df, created_gtf_df):
    pdt.assert_frame_equal(saved_gtf_df, created_gtf_df)


@pytest.fixture
def gtf_df_missing_features(created_gtf_df) -> pd.DataFrame:
    return create_missing_features(created_gtf_df, {})


def test_create_missing_features_identity(created_gtf_df, gtf_df_missing_features):
    pdt.assert_frame_equal(created_gtf_df, gtf_df_missing_features)

@pytest.fixture
def gtf_df_with_created_features(created_gtf_df) -> pd.DataFrame:
    return create_missing_features(
        created_gtf_df,
        unique_keys={"gene": "gene_id", "transcript": "transcript_id"},
        extra_columns={
            "gene": {"gene_name"},
            "transcript": {"gene_id", "gene_name", "transcript_name"},
        },
    )


def test_check_expanded_dataframe(gtf_df_with_created_features):
    df = gtf_df_with_created_features
    assert "gene" in df["feature"].values
    assert "transcript" in df["feature"].values

    C18orf10_201_transcript_mask = (df["feature"] == "transcript") & (  # noqa: N806
        df["transcript_name"] == "C18orf10-201"
    )
    C18orf10_201_df = df[C18orf10_201_transcript_mask]  # noqa: N806

    check_dataframe_values(df=C18orf10_201_df, seqname="18", start=32630766, end=32630768, strand="-")
    KIAA1328_df = df[(df["feature"] == "gene") & (df["gene_name"] == "KIAA1328")]  # noqa: N806
    check_dataframe_values(df=KIAA1328_df, seqname="18", start=32663078, end=32663157, strand="+")


def check_dataframe_values(df, seqname, start, end, strand):
    assert len(df) == 1
    assert df["seqname"].values[0] == seqname
    assert df["start"].values[0] == start
    assert df["end"].values[0] == end
    assert df["strand"].values[0] == strand


def test_missing_features(gtf_df_missing_features):
    assert "gene" not in gtf_df_missing_features["feature"].values
    assert "transcript" not in gtf_df_missing_features["feature"].values
