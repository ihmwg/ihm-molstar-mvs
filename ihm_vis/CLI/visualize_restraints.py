from ihm_vis import IHM_Builder
from ihm_vis import style
import argparse
import json

def visualize(cif_file, filters, sub_style_mode, max_restraints, style_file, output, title, verbose=False):

    # Set user style if provided
    if style_file is not None:
        if verbose: print(f"Setting style from: {style_file}")
        style.set_style_from_file(style_file)

    # Read in cif file from url or file
    ihm_b = IHM_Builder(cif_file)

    if ihm_b.source_type == "file" and verbose:
        print(f"Initialized IHM_Builder from local file")

    if ihm_b.source_type == "url" and verbose:
        print(f"Initialized IHM_Builder from url")

    # Parse all restraint types
    ihm_b.get_all_restraints()
    if verbose: print(f"Read {len(ihm_b.restraint_df)} restraints")

    # Apply user-requested filters
    for _filter in filters:
        len_before = len(ihm_b.restraint_df)

        if _filter == "random_sample" or _filter == "first_n":
            if verbose: print(f"Applying {_filter} filter with n={max_restraints}")
            ihm_b.filter_restraints(_filter, n=max_restraints)

        else:
            if verbose: print(f"Applying {_filter} filter")
            ihm_b.filter_restraints(_filter)

        len_after = len(ihm_b.restraint_df)
        if verbose:
            print(f"\tremoved {len_before - len_after} restraints")

    # Further sample to max restraints if necessary
    if len(ihm_b.restraint_df) > max_restraints:
        if verbose: print(f"Randomly sampling remaining restraints ({len(ihm_b.restraint_df)}) to requested max_restraints ({max_restraints})")
        ihm_b.filter_restraints("random_sample", n=max_restraints, random_state=27)

    # if sub_style mode is requested
    if sub_style_mode is not None:
        ihm_b.apply_sub_styles(sub_style_mode)
        if verbose: 
            print(f"Applied sub_style_mode {sub_style_mode} resulting in...")
            sub_style_outcomes = ihm_b.restraint_df["sub_style"].value_counts().to_dict()
            for ss, count in sub_style_outcomes.items():
                print(f"\t{ss}: {count} restraints")

    # Set restraint styles
    ihm_b.set_all_restraint_styles()
    if verbose: print(f"Applied styles to {len(ihm_b.restraint_df)} restraints")

    # Create output visualizatoin
    if output is None:
        output = ihm_b.basename.replace(".cif", ".mvsj")

    if title is None:
        title = ihm_b.basename

    if verbose: print(f"Writing output to {output} with title {title}")
    ihm_b.to_mvsj(output, title=title)

    if verbose: print("DONE")


def main():

    # Some filter funcs only make sense
    # to call from code because they require
    # additional argumets - such as the 'within_X_of' filter
    # only offer the user options that make sense for CLI
    # could add a lot of additional arguments to take in this input in future
    CLI_supported_filter_funcs = ["first_n", "random_sample", 
                                  "across_chains", "within_chains", 
                                  "violated", "compliant"]


    # Similar situation for sub_style_modes
    CLI_supported_sub_style_modes = ["violated_and_compliant", "inter_and_intra_chain"]

    ###########################################################################
    # Parse args
    ############

    parser = argparse.ArgumentParser(description="Visualize CIF file structures with IHM restraint data using ihm_vis.")
    parser.add_argument("cif_file",  
                        type=str, 
                        help="Path or URL of CIF file")

    parser.add_argument("-n", 
                        "--max_restraints", 
                        type=int, 
                        default=50, 
                        help="Maximum number of restraints to visualize, this value will be used if the 'first_n' or 'random_sample' filters are applied.")

    parser.add_argument("-f", 
                        "--filter", 
                        nargs="+",
                        action="append",
                        default=[],
                        choices=CLI_supported_filter_funcs, 
                        help=f"Filter(s) to apply to restraint table. You can specify multiple to chain filters together.")

    parser.add_argument("-s", 
                        "--sub_style_mode", 
                        choices=CLI_supported_sub_style_modes, 
                        default=None, 
                        help="Sub-style mode for coloring on a per-restraint basis")

    parser.add_argument("-c", 
                        "--style_file", 
                        type=str, 
                        default=None, 
                        help="Optional YAML or JSON file to customize visual style of macromolecule, restraint residues, or restraint distance primitives")

    parser.add_argument("-o", 
                        "--output", 
                        type=str,
                        default=None, 
                        help="Output file name. Note! when using a local file a .mvsx archive will be created bundling the .mvsj file with the CIF file. Defaults to CIF file stem.")

    parser.add_argument("-t", 
                        "--title", 
                        type=str, 
                        default=None, 
                        help="Title for mvsj file. Defaults to input name of input CIF file.")

    parser.add_argument("-v",
                       "--verbose",
                       action="store_true",
                       default=False,
                       help="Verbose outputs")


    args = parser.parse_args()
    all_filters = []
    for filter_sublist in args.filter:
        all_filters.extend(filter_sublist)


    ###########################################################################
    # Driver func
    #############

    visualize(
        cif_file=args.cif_file,
        filters=all_filters,
        sub_style_mode=args.sub_style_mode,
        max_restraints=args.max_restraints,
        style_file=args.style_file,
        output=args.output,
        title=args.title,
        verbose=args.verbose,
    )
