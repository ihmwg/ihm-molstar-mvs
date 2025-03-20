"""
Functions for parsing cif files and extracting restraint information
"""

from mmcif.io.PdbxReader import PdbxReader
import pandas as pd


def parse_file(filepath):
    """
    Parse an mmCIF (macromolecular Crystallographic Information File) file and return its contents.

    Parameters 
    ----------
    filepath : str
        The path to the mmCIF file that will be parsed.

    Returns
    -------
    containers : list
        A list containg the parsed data from the mmCIF file. 
    """
    # list to hold all data
    containers = []
    # open the mmcif file
    with open(filepath, 'r') as ifh:
        # load data into containers list
        PdbxReader(ifh).read(containers)
    return containers


def get_restraints(containers):
    """
    Extract cross-link restraint data from a list of containers and return it as a pandas DataFrame.

    Parameters
    ----------
    containers : list
        A list of containers, each containing objects with cross-link restraint data.

    Returns
    -------
    df : pd.DataFrame
        A pandas DataFrame where each row represents a cross-link restraint, ith columns including:
        'entity_id_1', 'asym_id_1', 'seq_id_1', 'comp_id_1', 'atom_id_1', 
        'entity_id_2', 'asym_id_2', 'seq_id_2', 'comp_id_2', 'atom_id_2', 
        'model_granularity', 'distance_threshold', and 'restraint_type'.
    """
    for container in containers:
        restraints = []
        # Extract cross-link restraint data from the container
        cross_link_data = container.getObj('ihm_cross_link_restraint')
        # Iterate through each row of cross-link data
        for i in range(cross_link_data.getRowCount()):
            data = {
            "entity_id_1":  cross_link_data.getValue("entity_id_1", i),
            "asym_id_1": cross_link_data.getValue("asym_id_1", i),
            "seq_id_1": cross_link_data.getValue("seq_id_1", i),
            "comp_id_1": cross_link_data.getValue("comp_id_1", i),
            "atom_id_1" : cross_link_data.getValue("atom_id_1", i),

            "entity_id_2": cross_link_data.getValue("entity_id_2", i),
            "asym_id_2": cross_link_data.getValue("asym_id_2", i),
            "seq_id_2": cross_link_data.getValue("seq_id_2", i),
            "comp_id_2": cross_link_data.getValue("comp_id_2", i),
            "atom_id_2" : cross_link_data.getValue("atom_id_2", i),

            "model_granularity": cross_link_data.getValue("model_granularity", i),
            "distance_threshold": cross_link_data.getValue("distance_threshold", i),
            "restraint_type": cross_link_data.getValue("restraint_type", i)
            }
            restraints.append(data)
    # Convert the list of restraints to a pandas DataFrame.
    df = pd.DataFrame(restraints)
    # If atom ids are not specified default to carbon alpha (CA)
    df['atom_id_1'] = df['atom_id_1'].str.replace('.', 'CA', regex=False)
    df['atom_id_2'] = df['atom_id_2'].str.replace('.', 'CA', regex=False)
    return df


