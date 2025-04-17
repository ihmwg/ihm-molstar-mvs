"""
"""

from typing import Dict, Callable
import pandas as pd

from ihm_vis.filters import get_solved_distance

###################################################################################################
BUILTIN_SUB_STYLE_FUNCS: Dict[str, Callable] = {}
###################################################################################################

def default(ihm_builder) -> pd.Series:
    return pd.Series("default", index=ihm_builder.restraint_df.index)

BUILTIN_SUB_STYLE_FUNCS["default"] = default

###################################################################################################

def violated_and_compliant(ihm_builder) -> pd.Series:
    """

    NOTE AND BE CLEAR THAT THIS DOES MODIFY IHM_Builder passed in!!

    """

    # Add compliant info if not already present
    if not "compliant" in ihm_builder.restraint_df.columns:
        ihm_builder.restraint_df = get_solved_distance(ihm_builder.restraint_df)

    # Sort values for consistency of residue colors
    # first coloring always seems to stick
    # and want violated to take precidence
    ihm_bulder.restraint_df = ihm_builder.restraint_df.sort_values("compliant")

    return ihm_builder.restraint_df["compliant"].map({True: "compliant", False: "violated"})

BUILTIN_SUB_STYLE_FUNCS["violated_and_compliant"] = violated_and_compliant



def inter_and_intra_chain(ihm_builder) -> pd.Series:
    """

    NOTE AND BE CLEAR THAT THIS DOES MODIFY IHM_Builder passed in!!

    """

    # Sort values for consistency of residue colors
    # first coloring always seems to stick
    # and want inter to take precidence
    ihm_builder.restraint_df["same_chain"] = ihm_builder.restraint_df["asym_id_1"] == ihm_builder.restraint_df["asym_id_2"]
    ihm_builder.restraint_df = ihm_builder.restraint_df.sort_values("same_chain")

    return ihm_builder.restraint_df["same_chain"].map({True: "intra", False: "inter"})

BUILTIN_SUB_STYLE_FUNCS["inter_and_intra_chain"] = inter_and_intra_chain

