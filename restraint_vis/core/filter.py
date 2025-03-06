"""
Tools for filtering restraint table
"""

from restraint_vis.utils.general import restraint_type_to_operator

import numpy as np
import pandas as pd

def first_n(df: pd.DataFrame, n: int=5) -> pd.DataFrame:
    """
    Keep first n restraints

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints
    n : int
        The number to keep

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe
    """
    return df.head(5)

def random_sample(df: pd.DataFrame, n: int=100, random_state: Optional[Int]=None) -> pd.DataFrame:
    """
    Randomly rample n restraints

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints
    n : int
        The number to sample
    random_state : int | None
        Random state passed to pd.DataFrame.sample

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe
    """
    return df.sample(min(len(df), n), random_state=random_state)


def across_chains(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to restraints accross chains

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe
    """
    return df.loc[df["asym_1"] != df["asym_2"]]


def diversity_filter(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter dataframe so that each unique residue with a restraint has at least one restraint shown, 
    avoiding more than one restraint per residue when possible

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe
    """

    df_a = df[["asym_1", "seq_id_1"]].rename({"asym_1":"asym", "seq_id_1":"seq_id"}, axis=1)
    df_a["original_index"] = df.index

    df_b = df[["asym_2", "seq_id_2"]].rename({"asym_2":"asym", "seq_id_2":"seq_id"}, axis=1)
    df_b["original_index"] = df.index
 
    stacked = pd.concat((df_a, df_b), ignore_index=True).sample(frac=1.0)
    to_keep = stacked.drop_duplicates(subset=["asym", "seq_id"])["original_index"].unique()

    return df.loc[df.index.isin(to_keep)]

def get_solved_distance(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute solved distance and whether restraint was compliant
    Adds two columns: 'solved_distance' (foat) and 'compliant' (bool)

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints

    Returns
    -------
    df : pd.DataFrame
        The updated dataframe
    """

    df["solved_distance"] = df.apply(lambda row: np.linalg.norm(np.array(row["atom_id_1_coords"]) - np.array(row["atom_id_2_coords"])), axis=1)
    df["compliant"] = df.apply(lambda row: restraint_type_to_operator(row["restraint_type"])(row["solved_distance"], row["distance_threshold"]), axis=1)

    return df

def violated(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to only violated restraints
    Adds two columns: 'solved_distance' (foat) and 'compliant' (bool)

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe with added columns
    """

    df = get_solved_distance(df)
    return df.loc[df["compliant"] == False]

def compliant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to only compliant restraints
    Adds two columns: 'solved_distance' (foat) and 'compliant' (bool)

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe with added columns
    """

    df = get_solved_distance(df)
    return df.loc[df["compliant"] == True]


def within_X_of(df: pd.DataFrame, X: float, position: np.ndarray | List) -> pd.DataFrame:

    """
    Filter to restraints where at least one residue is with X angstroms of the given position

    Parameters 
    ----------
    df : pd.DataFrame
        The input dataframe of restraints

    X : float
        The desired distance threshold

    position: np.ndarray | List

    Returns
    -------
    df : pd.DataFrame
        The filtered dataframe
    """

    position = np.array(position)
    if position.shape != (3,):
        raise ValueError("Position must be a 1D array of size 3")

    res1_dist = df.apply(lambda row: np.linalg.norm(position - np.array(row["atom_id_1_coords"])), axis=1)
    res2_dist = df.apply(lambda row: np.linalg.norm(position - np.array(row["atom_id_2_coords"])), axis=1)

    return df.loc[(res1_dist < X) | (res2_dist < X)]


