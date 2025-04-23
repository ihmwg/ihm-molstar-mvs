from ihm_vis import IHM_Builder
from ihm_vis import style
import argparse
import json

# cif_file = "https://pdb-ihm.org/cif/9a3v.cif"
# filters = ["across_chains"]
# sub_style_mode = "violated_and_compliant"
# max_restraints = 50
# style_yaml = "test_user_config.yaml"  # Sets residue in violated restraints to red!

output_file_stem = "9a3v_inter-chain_restraints"

def main(cif_file, filters, sub_style_mode, max_restraints, style_yaml, output_file_stem, config_json):
    if config_json is not None:
        with open(config_json, "r") as f:
            config = json.load(f)
        max_restraints = config.get("n", max_restraints)



    # Set user style if provided
    if style_yaml is not None:
        style.set_style_from_yaml(style_yaml)

    # Read in cif file from url or file
    ihm_b = IHM_Builder(cif_file)

    # Parse all restraint types
    ihm_b.get_all_restraints()

    # Apply user-requested filters
    for _filter in filters:
        ihm_b.filter_restraints(_filter)

    # Sample to max restraints
    ihm_b.filter_restraints("random_sample", n=max_restraints, random_state=27)

    # if sub_style mode is requested
    if sub_style_mode is not None:
        ihm_b.apply_sub_styles(sub_style_mode)

    # Set restraint styles
    ihm_b.set_all_restraint_styles()

    # Create output visualizatoin
    ihm_b.to_mvsj(output_file_stem, title="minimal example nb")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize filtered IHM restraints using ihm_vis.")

    parser.add_argument("cif_file", type=str, help="Path or URL to CIF file")
    parser.add_argument("--n", type=int, default=5, help="Maximum number of restraints to visualize (default: 50)")
    parser.add_argument("--filter", choices=["across_chains", "within_chain"], default="across_chains",
                        help="Filter type to apply to restraints (default: across_chains)")
    parser.add_argument("--sub_style", choices=["violated", "compliant", "violated_and_compliant", None], default="violated_and_compliant",
                        help="Sub-style mode for coloring based on compliance")
    parser.add_argument("--style_yaml", type=str, default=None, help="Optional YAML file to override visual style")
    parser.add_argument("--output", type=str, default="ihm_restraints", help="Output file stem (default: ihm_restraints)")
    parser.add_argument("--JSONconfig", type=str, default=None, help="Optional JSON config file with advanced settings")

    args = parser.parse_args()

    main(
        cif_file=args.cif_file,
        filters=[args.filter],
        sub_style_mode=args.sub_style,
        max_restraints=args.n,
        style_yaml=args.style_yaml,
        output_file_stem=args.output,
        config_json=args.JSONconfig,
    )
