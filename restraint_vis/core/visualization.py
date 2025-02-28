"""
Underlying functions for visualizing with molviewspec
"""

from typing import Dict

from restraint_vis.utils import color_utils
from restraint_vis.utils.general import RESTRAINT_TYPE_TO_SYMBOL, restraint_type_to_symbol

import molviewspec as mvs
from molviewspec.nodes import ComponentExpression
import pandas as pd


# will prob to refactor this later,
# cif parameter only used if color=chains
# so prob would want fix that later
def visualize_macromolecule(structure, # molviewspec scene builder
                            cif,       # cif file info

                            color: str="chains",
                            representation: str="cartoon",
                            ): 
    """
    Visualize the macromolecule

    Parameters
    ----------
    structure : mvs.Structure, The structure to visualize
    color: str | Dict[str, str], the color to make the macromolecule

        string options:
            'chains' - color each chain a unique color

            a hex color code - a valid hex color code

            CSS4 color names - any of the CSS4 color names from matplotlib
                such as 'red', 'grey', ''
    Raises
    ______
        ValueError: color not found
    """


    # Chains
    if color == "chains":
        chain_info = cif.getObj("struct_asym")
        id_col_index = chain_info.getAttributeIndexDict()["id"]
        chains = chain_info.getColumn(id_col_index)

        for chain, hex_color in zip(chains, color_utils.get_n_colors(len(chains))):
            structure.component(selector=ComponentExpression(label_asym_id=chain)).representation(type=representation).color(color=hex_color)


    # otherwise, try to get from matplotlib
    else:
        hex_color = color_utils.get_color(color)
        polymer = structure.component(selector="polymer").representation(type=representation).color(hex_color)


# TODO: feels clunky to use, think about improvements
# do we need some abstraction/ shorthand or are we fine
# with a large argument list
def visualize_restraint(structure,

                        start_asym_id: str | int, 
                        start_seq_id: int, 

                        end_asym_id: str | int,
                        end_seq_id: int, 

                        distance: float,
                        restraint_type: str, 

                        representation: str="ball_and_stick",
                        residue_color: str="red", 
                        start_atom_id: str="CA", end_atom_id: str="CA",

                        radius: float=0.1, 
                        dash_length: float=0.1,
                        line_color=None,

                        label_template="Solved Distance: {{{{distance}}}}, Restraint Distance: {restraint_symbol}{distance}",
                        label_color=None,

                        focus: bool=True,
                        ):

    f"""
    Visualize a restraint

    Parameters
    ----------

        start_asym_id: str | int, asym_id of the starting residue
        start_seq_id: int, seq_id of the starting residue

        end_asym_id: str | int, asym_id of the ending residue
        end_seq_id: int, seq_id of the ending residue

        distance: float, the restraint distance
        restraint_type: str, the type/operator of the restraint, one of {list(RESTRAINT_TYPE_TO_SYMBOL.keys())}

        representation: str="ball_and_stick", representation of the two residues of the restraint
        residue_color: str="red", color of the restraint residues
        start_atom_id: str="CA", atom_id of the starting residue
        end_atom_id: str="CA", atom_id of the ending residue

        radius: float=0.1, radius of the line drawn between restraint residues
        line_color=None, color for the line drawn between restraint residues, defaults to the residues' color
        dash_length: float=0.1, dash length of the line drawn between restraint residues

        label_template="Solved Distance: {{{{distance}}}}, Restraint Distance: the text displayed on the line drawn between restraint residues
        label_color=None, color of the text displayed on the line drawn between residues, defaults to the residues' color

        focus: bool=True, whether to focus the carmera to this restraint
 
    Raises
    ______
        ValueError: color not found
    """


    # TODO: Cynthia's parser already handles this
    # decide if we should check twice or not
    if start_atom_id == ".":
        start_atom_id = "CA"

    if end_atom_id == ".":
        end_atom_id = "CA"

    residue_hex_color = color_utils.get_color(residue_color)

    if label_color is None:
        label_hex_color = residue_hex_color
    else:
        label_hex_color = color_utils.get_color(label_color)

    if line_color is None:
        line_hex_color = residue_hex_color
    else:
        line_hex_color = color_utils.get_color(line_color)


    start_residue = ComponentExpression(label_asym_id=start_asym_id,
                                       beg_label_seq_id=start_seq_id,
                                       end_label_seq_id=start_seq_id)

    end_residue = ComponentExpression(label_asym_id=end_asym_id,
                                       beg_label_seq_id=end_seq_id,
                                       end_label_seq_id=end_seq_id)

    start_atom = ComponentExpression(label_asym_id=start_asym_id,
                                       beg_label_seq_id=start_seq_id,
                                       end_label_seq_id=start_seq_id,
                                       label_atom_id=start_atom_id)

    end_atom = ComponentExpression(label_asym_id=end_asym_id,
                                       beg_label_seq_id=end_seq_id,
                                       end_label_seq_id=end_seq_id,
                                       label_atom_id=end_atom_id)



    start_component = structure.component(selector=start_residue)
    start_component.representation(type=representation).color(color=residue_hex_color)

    end_component = structure.component(selector=end_residue)
    end_component.representation(type=representation).color(color=residue_hex_color)

    res = structure.primitives().distance(
            start=start_atom,
            end=end_atom,
            color=line_hex_color,
            radius=radius,
            dash_length=dash_length,
            label_template=label_template.format(restraint_symbol=restraint_type_to_symbol(restraint_type), distance=distance),
            label_color=label_hex_color)

    if focus:
        res.focus()


