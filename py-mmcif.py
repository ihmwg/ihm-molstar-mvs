import argparse
from mmcif.io.PdbxReader import PdbxReader

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
        cross_link_data = container.getObj('ihm_cross_link_restraint')
        for i in range(cross_link_data.getRowCount()):
            entity_id_1 = cross_link_data.getValue("entity_id_1", i)
            asym_id_1 = cross_link_data.getValue("asym_id_1", i)
            seq_id_1 = cross_link_data.getValue("seq_id_1", i)
            comp_id_1 = cross_link_data.getValue("comp_id_1", i)

            entity_id_2 = cross_link_data.getValue("entity_id_2", i)
            asym_id_2 = cross_link_data.getValue("asym_id_2", i)
            seq_id_2 = cross_link_data.getValue("seq_id_2", i)
            comp_id_2 = cross_link_data.getValue("comp_id_2", i)

            model_granularity = cross_link_data.getValue("model_granularity", i)
            distance_threshold = cross_link_data.getValue("distance_threshold", i)
            restraint_type = cross_link_data.getValue("restraint_type", i)

            # Print extracted data
            print(f"Cross-link Restraint Data:")
            print(f"Entity 1: {entity_id_1}, Asym ID 1: {asym_id_1}, Seq ID 1: {seq_id_1}, Comp ID 1: {comp_id_1}")
            print(f"Entity 2: {entity_id_2}, Asym ID 2: {asym_id_2}, Seq ID 2: {seq_id_2}, Comp ID 2: {comp_id_2}")
            print(f"Model Granularity: {model_granularity}, Distance Threshold: {distance_threshold}, Restraint Type: {restraint_type}")
            print("-" * 50)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filepath", help='Path to mmcif file')
    args = parser.parse_args()

    containers = parse_file(args.filepath)
    get_restraints(containers)
    