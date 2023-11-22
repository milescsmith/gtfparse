from importlib.resources import as_file, files

import pandas as pd
import pytest

from gtfparse.read_gtf import read_gtf

# ruff: noqa: S101


@pytest.fixture
def b16_string_gtf():
    with as_file(files("tests.data").joinpath("B16.stringtie.head.gtf")) as gtf:
        return read_gtf(gtf, expand_attribute_column=True)


def check_b16_required_columns(b16_string_gtf):
    assert pd.Series(["feature", "cov", "FPKM"]).isin(b16_string_gtf.columns).all()


def check_b16_features(b16_string_gtf):
    assert pd.Series(["exon", "transcript"]).isin(b16_string_gtf["feature"]).all()


def check_string_cov_and_FPKM(b16_string_gtf):  # noqa: N802
    for i in b16_string_gtf.itertuples():
        if i.feature == "exon":
            assert float(i.cov) >= 0
        elif i.feature == "transcript":
            assert float(i.cov) >= 0
            assert float(i.FPKM) >= 0


@pytest.fixture
def b16_float_gtf():
    with as_file(files("tests.data").joinpath("B16.stringtie.head.gtf")) as gtf:
        return read_gtf(gtf, expand_attribute_column=True, column_converters={"cov": float, "FPKM": float})


def check_float_cov_and_FPKM(b16_float_gtf):  # noqa: N802
    for i in b16_float_gtf.itertuples():
        assert isinstance(i.cov, float)
        if i.feature == "exon":
            assert i.cov >= 0
        elif i.feature == "transcript":
            assert isinstance(i.FPKM, float)
            assert i.cov >= 0
            assert i.FPKM >= 0
