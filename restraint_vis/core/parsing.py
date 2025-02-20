"""
Functions for parsing cif files and extracting restraint information
"""

from mmcif.io.PdbxReader import PdbxReader
import pandas as pd


def parse_file(filepath):
    # list to hold all data
    containers = []
    # open the mmcif file
    with open(filepath, 'r') as ifh:
        # load data into containers list
        PdbxReader(ifh).read(containers)
    return containers


def get_restraints(containers):
    for container in containers:
        restraints = []
        cross_link_data = container.getObj('ihm_cross_link_restraint')
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

    df = pd.DataFrame(restraints)
    df['atom_id_1'] = df['atom_id_1'].str.replace('.', 'CA', regex=False)
    df['atom_id_2'] = df['atom_id_2'].str.replace('.', 'CA', regex=False)
    return df


