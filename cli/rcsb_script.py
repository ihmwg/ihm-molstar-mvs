import molviewspec as mvs
from molviewspec.nodes import ComponentExpression
import mmcif
from mmcif.io.PdbxReader import PdbxReader
import pandas as pd
import argparse
import os
import sys
from restraint_vis import core

def parse_file(filepath):
    containers = []
    with open(filepath, 'r') as ifh:
        PdbxReader(ifh).read(containers)
        return containers[0]


def get_restraints(cif):
    restraints = []
    cross_link_data = cif.getObj('ihm_cross_link_restraint')
    
    for i in range(cross_link_data.getRowCount()):
        data = {
                "asym_id_1": cross_link_data.getValue("asym_id_1", i),
                "seq_id_1": cross_link_data.getValue("seq_id_1", i),
                "atom_id_1" : cross_link_data.getValue("atom_id_1", i),
    
                "asym_id_2": cross_link_data.getValue("asym_id_2", i),
                "seq_id_2": cross_link_data.getValue("seq_id_2", i),
                "atom_id_2" : cross_link_data.getValue("atom_id_2", i),
            
                "distance_threshold": cross_link_data.getValue("distance_threshold", i),
                "restraint_type": cross_link_data.getValue("restraint_type", i),
            }
    
        restraints.append(data)
        df = pd.DataFrame(restraints)
        df['atom_id_1'] = df['atom_id_1'].str.replace('.', 'CA', regex=False)
        df['atom_id_2'] = df['atom_id_2'].str.replace('.', 'CA', regex=False)
        
    return df

def get_mvsj(filepath, parsed_cif, restraints_df):
    builder = mvs.create_builder()
    structure = builder.download(url=f"https://pdb-ihm.org/cif/{filepath}").parse(format="mmcif").assembly_structure(assembly_id="1")
    core.visualize_macromolecule(structure, parsed_cif)

    row = restraints_df.iloc[0]

    restraint_info = dict(
        start_asym_id = row["asym_id_1"],
        start_seq_id = row["seq_id_1"],
        start_atom_id = row["atom_id_1"],
        
        end_asym_id = row["asym_id_2"],
        end_seq_id = row["seq_id_2"],
        end_atom_id = row["atom_id_2"],
        
        distance = row["distance_threshold"],
        restraint_type = row["restraint_type"],
    )

    core.visualize_restraint(structure, **restraint_info)

    builder.save_state(destination=f"{filepath[:4]}_restraint.mvsj", title="with restraint")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Pass a cif file through the command-line.")

    parser.add_argument("cif_file_name", type=str, help="The CIF file to process")

    args = parser.parse_args()

    parsed_cif = parse_file(args.cif_file_name)

    restraint_pd = get_restraints(parsed_cif)

    get_mvsj(args.cif_file_name, parsed_cif, restraint_pd)
