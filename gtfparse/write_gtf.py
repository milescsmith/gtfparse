import pandas as pd
from .required_columns import REQUIRED_COLUMNS


def extract_seq_info(gtf_row: pd.Series):
    """ Convert the data in the rows from a GTF/GFF3-formatted DataFrame
    into the appropriate tab- and semicolon-delimited string representations
    \f
    Parameters
    ----------

    gtf_row : :class:`pd.Series`
        Row of a GTF-formatted Pandas DataFrame
    """
    try:
        required = "\t".join(list(gtf_row[REQUIRED_COLUMNS[:-1]].astype("str").values))
        attributes = ";".join(
            [
                "=".join(_)
                for _ in zip(
                    gtf_row.drop(REQUIRED_COLUMNS[:-1]).keys().values,
                    list(gtf_row.drop(REQUIRED_COLUMNS[:-1]).astype("str")),
                )
            ]
        )
        line_to_strings = required + "\t" + attributes + "\n"
    except Exception:
        raise ValueError(gtf_row)
    return line_to_strings


def df_to_gtf(df: pd.DataFrame, filename: str) -> None:
    """ Write a GTF/GFF3-formatted DataFrame out as a GTF
    \f
    Parameters
    ----------

    df : :class:`pd.DataFrame`
        GTF in the form of a Pandas DataFrame

    filename : `str`
        name of file to write GTF out as
    """
    line_to_strings = df.apply(extract_seq_info, axis=1)
    with open(filename, "w") as gtfoutput:
        # gtfoutput.writelines("seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n")
        gtfoutput.writelines(line_to_strings)
