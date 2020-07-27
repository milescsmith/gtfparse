# Copyright (c) 2015-2018. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
from math import ceil
from os import stat
from os.path import exists, basename
from sys import intern
from typing import Optional, List, Callable, Dict, Union

import numpy as np
import pandas as pd
from tqdm import tqdm

from .attribute_parsing import expand_attribute_strings
from .parsing_error import ParsingError
from .logging import setup_logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_gtf(
    filepath_or_buffer: str,
    chunksize: int = 1024 * 1024,
    features: Optional[List[str]] = None,
    intern_columns: List[str] = ["seqname", "source", "strand", "frame"],
    fix_quotes_columns: List[str] = ["attribute"],
):
    """
    Parameters
    ----------

    filepath_or_buffer : str or buffer object

    chunksize : int

    features : set or None
        Drop entries which aren't one of these features

    intern_columns : list
        These columns are short strings which should be interned

    fix_quotes_columns : list
        Most commonly the 'attribute' column which had broken quotes on
        some Ensembl release GTF files.
    """
    logger = logging.getLogger("gtfparse")

    if features:
        features = np.unique(features)

    dataframes = []

    # GTF columns:
    # 1) seqname: str ("1", "X", "chrX", etc...)
    # 2) source : str
    #      Different versions of GTF use second column as of:
    #      (a) gene biotype
    #      (b) transcript biotype
    #      (c) the annotation source
    #      See: https://www.biostars.org/p/120306/#120321
    # 3) feature : str ("gene", "transcript", &c)
    # 4) start : int
    # 5) end : int
    # 6) score : float or "."
    # 7) strand : "+", "-", or "."
    # 8) frame : 0, 1, 2 or "."
    # 9) attribute : key-value pairs separated by semicolons
    # (see more complete description in docstring at top of file)

    def parse_frame(s: str) -> int:
        if s == ".":
            return 0
        else:
            return int(s)

    def fix_attribute_column(x: str) -> str:
        return x.replace(';"', '"').replace(";-", "-").replace("; ", ";")

    # tqdm.pandas(tqdm, leave=True)
    logger.info(f"Reading in {basename(filepath_or_buffer)} chunks")

    chunk_iterator = pd.read_csv(
        filepath_or_buffer,
        sep="\t",
        comment="#",
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
        dtype={"score": np.float32, "attribue": str, "strand": "category"},
        skipinitialspace=True,
        skip_blank_lines=True,
        error_bad_lines=True,
        warn_bad_lines=True,
        chunksize=chunksize,
        engine="c",
        na_values=".",
        converters={
            "frame": parse_frame,
            "seqname": intern,
            "source": intern,
            "feature": intern,
        },
        memory_map=True,
        low_memory=False,
    )

    file_size = stat(filepath_or_buffer).st_size

    try:
        for df in tqdm(
            chunk_iterator,
            desc="loading file",
            total=ceil(file_size / (chunksize * 425)),
            unit="chunks",
            leave=True,
        ):
            dataframes.append(df)
    except Exception as e:
        raise ParsingError(str(e))

    df = pd.concat(dataframes)
    if features:
        logger.info(f"Filtering for entries that have a feature in {features}")
        df = df[df["feature"].isin(features)]

    try:
        import swifter

        logger.info("swifter found, processing in parallel")
        logger.info("Repairing non-standard 'attributes'")
        df["attribute"] = (
            df["attribute"]
            .swifter.progress_bar(True)
            .apply(lambda x: x.replace(';"', '"').replace(";-", "-").replace("; ", ";"))
        )

        logger.info("Converting non-integer 'start' values to 0")
        df["start"] = (
            df["start"]
            .swifter.progress_bar(True)
            .apply(lambda i: np.nan_to_num(i))
            .astype(np.int32)
        )
        logger.info("Converting non-integer 'end' values to 0")
        df["end"] = (
            df["end"]
            .swifter.progress_bar(True)
            .apply(lambda i: np.nan_to_num(i))
            .astype(np.int32)
        )
    except ImportError:
        logger.info("Repairing non-standard 'attributes'")
        df["attribute"] = df["attribute"].apply(
            lambda x: x.replace(';"', '"').replace(";-", "-").replace("; ", ";")
        )
        logger.info("Converting non-integer 'start' values to 0")
        df["start"] = df["start"].apply(lambda i: np.nan_to_num(i)).astype(np.int32)
        logger.info("Converting non-integer 'end' values to 0")
        df["end"] = df["end"].apply(lambda i: np.nan_to_num(i)).astype(np.int32)

    return df


def parse_gtf_and_expand_attributes(
    filepath_or_buffer: str,
    chunksize: int = 1024 * 1024,
    restrict_attribute_columns: Union[List[str]] = None,
    features: List[str] = None,
):
    """
    Parse lines into column->values dictionary and then expand
    the 'attribute' column into multiple columns. This expansion happens
    by replacing strings of semi-colon separated key-value values in the
    'attribute' column with one column per distinct key, with a list of
    values for each row (using None for rows where key didn't occur).

    Parameters
    ----------
    filepath_or_buffer : str or buffer object

    chunksize : int

    restrict_attribute_columns : list/set of str or None
        If given, then only attribute columns.

    features : set or None
        Ignore entries which don't correspond to one of the supplied features
    """
    logger = logging.getLogger("gtfparse")

    result = parse_gtf(filepath_or_buffer, chunksize=chunksize, features=features)

    logger.info("Expanding attributes")
    expanded_result = expand_attributes(
        df=result, restrict_attribute_columns=restrict_attribute_columns
    )

    return expanded_result


def expand_attributes(
    df: pd.DataFrame, restrict_attribute_columns: Optional[List[str]] = None
) -> pd.DataFrame:

    logger = logging.getLogger("gtfparse")

    def attribute_to_dict(
        x: str, delim1: str = ";", delim2: str = "="
    ) -> Dict[str, str]:
        return dict(
            _.split(delim2, maxsplit=1)
            for _ in x.split(delim1)
            if len(_.split(delim2, maxsplit=1)) == 2
        )

    try:
        import swifter

        logger.info("Converting 'attribute' column to dictionaries, then to json")
        attribute_values = pd.read_json(
            df["attribute"]
            .swifter.progress_bar(True)
            .apply(attribute_to_dict)
            .to_json(orient="records")
        )
    except ImportError:
        logger.info("Converting 'attribute' column to dictionaries, then to json")
        attribute_values = pd.read_json(
            df["attribute"].apply(attribute_to_dict).to_json(orient="records")
        )

    if restrict_attribute_columns:
        logger.info("Combining 'attribute' values not selected for expansion")
        fold_columns = attribute_values.columns[
            ~attribute_values.columns.isin(restrict_attribute_columns)
        ]
        attribute_column = attribute_values[fold_columns].apply(
            lambda x: ";".join(f"{k}={v}" for k, v in x.items() if not pd.isnull(v)),
            axis=1,
        )
        logger.info("Concatenating columns")
        expanded_df = pd.concat(
            [
                df.loc[:, df.columns.drop("attribute")],
                attribute_column,
                attribute_values,
            ],
            axis=1,
        )
    else:
        logger.info("Concatenating columns")
        expanded_df = pd.concat(
            [df.loc[:, df.columns.drop("attribute")], attribute_values], axis=1
        )

    return expanded_df


def read_gtf(
    filepath_or_buffer: str,
    expand_attribute_column: bool = False,
    infer_biotype_column: bool = False,
    column_converters: Optional[Dict[str, Callable[..., str]]] = None,
    usecols: Optional[List[str]] = None,
    features: List[str] = None,
    chunksize: int = 1024 * 1024,
):
    """
    Parse a GTF into a dictionary mapping column names to sequences of values.

    Parameters
    ----------
    filepath_or_buffer : str or buffer object
        Path to GTF file (may be gzip compressed) or buffer object
        such as StringIO

    expand_attribute_column : bool
        Replace strings of semi-colon separated key-value values in the
        'attribute' column with one column per distinct key, with a list of
        values for each row (using None for rows where key didn't occur).

    infer_biotype_column : bool
        Due to the annoying ambiguity of the second GTF column across multiple
        Ensembl releases, figure out if an older GTF's source column is actually
        the gene_biotype or transcript_biotype.

    column_converters : dict, optional
        Dictionary mapping column names to conversion functions. Will replace
        empty strings with None and otherwise passes them to given conversion
        function.

    usecols : list of str or None
        Restrict which columns are loaded to the give set. If None, then
        load all columns.

    features : set of str or None
        Drop rows which aren't one of the features in the supplied set

    chunksize : int
    """
    setup_logging(name="gtfparse")
    
    logger = logging.getLogger("gtfparse")

    if isinstance(filepath_or_buffer, str) and not exists(filepath_or_buffer):
        logger.exception(f"GTF file does not exist: {filepath_or_buffer}")
        raise ValueError

    if expand_attribute_column:
        result_df = parse_gtf_and_expand_attributes(
            filepath_or_buffer, chunksize=chunksize, restrict_attribute_columns=usecols
        )
    else:
        result_df = parse_gtf(
            filepath_or_buffer, chunksize=chunksize, features=features
        )

    if column_converters:
        for column_name, column_type in list(column_converters.items()):
            result_df[column_name] = [
                column_type(string_value) if string_value else None
                for string_value in result_df[column_name]
            ]

    # Hackishly infer whether the values in the 'source' column of this GTF
    # are actually representing a biotype by checking for the most common
    # gene_biotype and transcript_biotype value 'protein_coding'
    if infer_biotype_column:
        unique_source_values = result_df["source"].unique()
        if "protein_coding" in unique_source_values:
            column_names = result_df.columns.unique()
            # Disambiguate between the two biotypes by checking if
            # gene_biotype is already present in another column. If it is,
            # the 2nd column is the transcript_biotype (otherwise, it's the
            # gene_biotype)
            if "gene_biotype" not in column_names:
                logger.info("Using column 'source' to replace missing 'gene_biotype'")
                result_df["gene_biotype"] = result_df["source"]
            if "transcript_biotype" not in column_names:
                logger.info(
                    "Using column 'source' to replace missing 'transcript_biotype'"
                )
                result_df["transcript_biotype"] = result_df["source"]

    if usecols is not None:
        column_names = result_df.columns.unique()
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    return result_df
