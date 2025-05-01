"""
"""

from typing import Dict, Callable
import pandas as pd

from ihm_vis.utils.restraint_filters import get_solved_distance

###################################################################################################
BUILTIN_SUB_STYLE_MODES: Dict[str, Callable] = {}
###################################################################################################

def default(df):
    df["sub_style"] = "default"
    return df

BUILTIN_SUB_STYLE_MODES["default"] = default

###################################################################################################

def violated_and_compliant(df):
    """
    """

    # Add compliant info if not already present
    if not "compliant" in df.columns:
        df = get_solved_distance(df)

    # Sort values for consistency of residue colors
    # first coloring always seems to stick
    # and want violated to take precidence
    df = df.sort_values("compliant", ascending=False)
    df["sub_style"] = df["compliant"].map({True: "compliant", False: "violated"})

    return df

BUILTIN_SUB_STYLE_MODES["violated_and_compliant"] = violated_and_compliant



def inter_and_intra_chain(df):
    """
    """

    # Sort values for consistency of residue colors
    # first coloring always seems to stick
    # and want inter to take precidence
    df["same_chain"] = df["asym_id_1"] == df["asym_id_2"]
    df = df.sort_values("same_chain", ascending=False)
    df["sub_style"] = df["same_chain"].map({True: "intra", False: "inter"})

    return df

BUILTIN_SUB_STYLE_MODES["inter_and_intra_chain"] = inter_and_intra_chain

