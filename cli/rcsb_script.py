import molviewspec as mvs
from molviewspec.nodes import ComponentExpression
import mmcif
from mmcif.io.PdbxReader import PdbxReader
import pandas as pd
import zipfile
import argparse
import sys
import os

def retrieve_cross_links_restraint_info(cif_file_path):
    if os.path.exists(cif_file_path):
        print(f"Found file: {cif_file_path}")
    else:
        print("File does not exist")
        exit(1)

    containers = []
    with open(cif_file_path, 'r') as f:
        mmcif.io.PdbxReader.PdbxReader(f).read(containers)

    cif = containers[0]

    cross_link_lists = cif.getObj('ihm_cross_link_list')
    cross_link_restraints = cif.getObj('ihm_cross_link_restraint')


    restraints = []

    for i in range(cross_link_restraints.getRowCount()):
        cross_link_data = {
                "id": cross_link_restraints.getValue("id", i),
                "entity_id_1": cross_link_restraints.getValue("entity_id_1", i),
                "entity_id_2": cross_link_restraints.getValue("entity_id_1", i),
                "asym_id_1": cross_link_restraints.getValue("asym_id_1", i),
                "asym_id_2": cross_link_restraints.getValue("asym_id_2", i),
                "seq_id_1": cross_link_restraints.getValue("seq_id_1", i),
                "comp_id_1": cross_link_restraints.getValue("comp_id_1", i),
                "seq_id_2": cross_link_restraints.getValue("seq_id_2", i),
                "comp_id_2": cross_link_restraints.getValue("comp_id_2", i),
                "atom_id_1" : cross_link_restraints.getValue("atom_id_1", i),
                "atom_id_2" : cross_link_restraints.getValue("atom_id_2", i),
                "distance_threshold": cross_link_restraints.getValue("distance_threshold", i),
                "restraint_type": cross_link_restraints.getValue("restraint_type", i),

            }

        restraints.append(cross_link_data)

    restraints_pd = pd.DataFrame(restraints)

    return restraints_pd

    # print(restraints)


def build_MVS_representation(cif_file_path, restraints):
    builder = mvs.create_builder()
    structure = builder.download(url=f"https://pdb-ihm.org/cif/{cif_file_path}").parse(format="mmcif").assembly_structure(assembly_id="1")
    polymer = structure.component(selector="polymer")
    representation = polymer.representation(type="cartoon")
    representation.color(color="#4b7fcc")

    restraint_type_to_repr = {
    "upper bound": "<",
    }

    for i, row in restraints.iterrows():
        start_entity_id = int(row["entity_id_1"])
        end_entity_id = int(row["entity_id_2"])
        start_seq_id = int(row["seq_id_1"])
        end_seq_id = int(row["seq_id_2"])
        start_atom_id = row["atom_id_1"]
        end_atom_id = row["atom_id_2"]
        threshold = float(row["distance_threshold"])
        restraint_type = row["restraint_type"]

        start_residue = ComponentExpression(label_entity_id=start_entity_id, 
                                            beg_label_seq_id=start_seq_id, end_label_seq_id=start_seq_id
                                        )
        start_atom = ComponentExpression(label_entity_id=start_entity_id, 
                                        beg_label_seq_id=start_seq_id, end_label_seq_id=start_seq_id,
                                        label_atom_id=start_atom_id
                                        )

        
        end_residue = ComponentExpression(label_entity_id=end_entity_id, 
                                        beg_label_seq_id=end_seq_id, end_label_seq_id=end_seq_id
                                        )

        end_atom = ComponentExpression(label_entity_id=end_entity_id, 
                                    beg_label_seq_id=end_seq_id, end_label_seq_id=end_seq_id, 
                                    label_atom_id=start_atom_id
                                    )

        structure.primitives().distance(
                start=start_atom,
                end=end_atom,
                color="red",
                radius=0.1,
                dash_length=0.1,
                label_template="Solved Distance: {{distance}}" + f", Restraint Distance: {restraint_type_to_repr[restraint_type]}{threshold}",
                label_color="red").focus()

        
        start_comp = structure.component(selector=start_residue)
        start_comp.representation(type="ball_and_stick").color(color="#ff0000")

        end_comp = structure.component(selector=end_residue)
        end_comp.representation(type="ball_and_stick").color(color="#ff0000")

    builder.save_state(destination="test.mvsj", title="8zz1_visualized", indent=2)

    # return builder.get_state()
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pass a cif file through the command-line.")

    parser.add_argument("cif_file_name", type=str, help="The CIF file to process")

    args = parser.parse_args()

    restraints_pd = retrieve_cross_links_restraint_info(args.cif_file_name)

    build_MVS_representation(args.cif_file_name, restraints_pd)

    # print(json_file)