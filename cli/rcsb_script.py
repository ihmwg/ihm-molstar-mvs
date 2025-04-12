import molviewspec as mvs
import pandas as pd
import argparse
import numpy as np
from restraint_vis import core, config
from restraint_vis.utils.general import restraint_type_to_operator
from pathlib import Path
import json

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
    return df.loc[df["asym_id_1"] != df["asym_id_2"]]

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
    return df.head(n)

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

def violated_and_compliant_mode(structure, df):
    df = core.filter_funcs.get_solved_distance(df)
    df["sub_style"] = df["compliant"].map({True: "compliant", False: "violated"})

    # Sort so that violated restraints are visualized last, 
    # so that if a residue is part of a compliant AND violated restraint
    # it is colored as violated
    df = df.sort_values("sub_style", ascending=False)
    # display(df)
    core.visualize_restraints(structure, df, sub_style_col="sub_style")


def inter_and_intra_chain_mode(structure, df):
    df["sub_style"] = (df["asym_id_1"] == df["asym_id_2"]).map({True: "intra", False: "inter"})

    # Sort so that inter restraints are visualized last, 
    # so that if a residue is part of a inter AND intra restraint
    # it is colored as inter
    df = df.sort_values("sub_style", ascending=False)
    core.visualize_restraints(structure, df, sub_style_col="sub_style")


def visualize_restraints(structure, restraint_df, sub_style_col=None, **kwargs):

    if sub_style_col is not None and not sub_style_col in restraint_df.columns:
        raise ValueError(f"The requested sub_style col ({sub_style_col}) is not present in the dataframe!")

    for _, row in restraint_df.iterrows():

        restraint_info = {
            "start_asym_id"  : row["asym_id_1"],
            "start_seq_id"   : row["seq_id_1"],
            "start_atom_id"  : row["atom_id_1"],

            "end_asym_id"    : row["asym_id_2"],
            "end_seq_id"     : row["seq_id_2"],
            "end_atom_id"    : row["atom_id_2"],

            "distance"       : row["distance_threshold"],
            "restraint_type" : row["restraint_type"],
        }

        if sub_style_col is not None:
            ss_val = row[sub_style_col]
        else:
            ss_val = None

        core.visualize_restraint(structure, **restraint_info, sub_style=ss_val, **kwargs)


def validate_and_prepare_inputs(params):

    cif_file = params["cif_file_name"]
    filters = params.get("filter", ["violated"])  # default fallback
    sub_style_mode = params.get("sub_style_mode")
    style_json = params.get("style_json")
    max_restraints = params.get("n", 10)
    output_file = params.get("output_file", f"{cif_file}_output.mvsj")

    for f in filters:
        if f not in FILTERS:
            raise ValueError(f"Unknown filter '{f}'. Valid options: {list(FILTERS.keys())}")
    
    if not Path(cif_file).exists():
        raise FileNotFoundError(f"The specified CIF file ({cif_file}) does not exist!")

    if sub_style_mode and sub_style_mode not in SUB_STYLE_MODES:
        raise ValueError(f"Unknown sub_style_mode '{sub_style_mode}'. Valid options: {list(SUB_STYLE_MODES.keys())}")

    if style_json:
        if not Path(style_json).exists():
            raise FileNotFoundError(f"Style JSON file '{style_json}' not found")
        config.set_style_from_json(style_json)

    return cif_file, filters, sub_style_mode, style_json, max_restraints, output_file

def get_mvsj(filepath, n, filter):

    cif = core.parse_file(filepath)
    builder = mvs.create_builder()
    structure = builder.download(url=f"https://pdb-ihm.org/cif/{filepath}").parse(format="mmcif").assembly_structure(assembly_id="1")
    core.visualize_macromolecule(structure, cif[0])

    # Parse and visualize restraints
    restraints = core.get_restraints(cif)
    restraints = get_solved_distance(restraints)
    # restraints = across_chains(restraints)

    filter = filter.lower()

    if filter == "compliant":
        restraints = compliant(restraints)
    elif filter == "violated":
        restraints = violated(restraints)

    restraints = first_n(restraints, n)
    print(restraints)

    visualize_restraints(structure, restraints)

    builder.save_state(destination=f"{filepath[:4]}_restraint.mvsj", title="with restraint")


def load_config(json_file):
    """ Load parameters from a JSON file. """
    with open(json_file, 'r') as f:
        return json.load(f)
    
def merge_args_and_config(args, config):
    merged = config.copy()
    for k, v in vars(args).items():
        if v is not None:
            merged[k] = v
    return merged
    

def run_pipeline(cif_file, filters, sub_style_mode, max_restraints, output_file):
    cif = core.parse_file(cif_file)
    restraint_df = core.get_restraints(cif)

    print("Restraints before filtering:")
    print(restraint_df)

    builder = mvs.create_builder()
    structure = builder.download(url=f"https://pdb-ihm.org/cif/{cif_file}").parse(format="mmcif").assembly_structure()
    core.visualize_macromolecule(structure, cif[0])

    for f in filters:
        n_before = len(restraint_df)
        restraint_df = FILTERS[f](restraint_df)
        print(f"{f} removed {n_before - len(restraint_df)} restraints")

    if len(restraint_df) > max_restraints:
        restraint_df = core.filter_funcs.random_sample(restraint_df, n=max_restraints)

    print("Restraints after filtering:")
    print(restraint_df)

    if sub_style_mode:
        SUB_STYLE_MODES[sub_style_mode](structure, restraint_df)
    else:
        core.visualize_restraints(structure, restraint_df)

    print(f"Saving MVSJ output to: {output_file}")
    builder.save_state(destination=output_file, title="with restraint")

    

    # Validate inputs
FILTERS = {
               "across_chains": across_chains,
               "diversity_filter": core.filter_funcs.diversity_filter,
               "violated": violated,
               "compliant": compliant,
           }
    
SUB_STYLE_MODES = {

            "violated_and_compliant": violated_and_compliant_mode,
            "inter_and_intra_chain": inter_and_intra_chain_mode,
    }


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Pass a cif file & number of restraints through the command-line.")

    parser.add_argument("cif_file_name", type=str, help="The CIF file to process")

    parser.add_argument("--n", type=int, default=5, help="Number of restraints to keep (default: 5)")

    parser.add_argument("--filter", nargs='+', choices=["violated", "compliant", "across_chains", "diversity_filter"],
                    default=["violated"])
    
    parser.add_argument("--chain", type=str, choices=["inter", "intra"], default="inter",
                        help="Filter restraints: 'Inter' or 'Intra' chain restraints" )
    
    parser.add_argument("--config", type=str, help="Path to a JSON config file")

    args = parser.parse_args()
    config_data = load_config(args.config) if args.config else {}
    merged_config = merge_args_and_config(args, config_data)

    cif_file, filters, sub_style_mode, style_json, max_restraints, output_file = validate_and_prepare_inputs(merged_config)

    print(f"\nFinal Configuration:\n"
          f"- CIF file: {cif_file}\n"
          f"- Filters: {filters}\n"
          f"- Max restraints: {max_restraints}\n"
          f"- Sub-style mode: {sub_style_mode}\n"
          f"- Style config: {style_json}\n"
          f"- Output file: {output_file}\n")

    run_pipeline(cif_file, filters, sub_style_mode, max_restraints, output_file)
