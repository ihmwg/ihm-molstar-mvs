"""
Underlying functions for visualizing with molviewspec
"""

from typing import Dict

from restraint_vis.utils import color_utils
from restraint_vis.utils.general import restraint_type_to_symbol

import molviewspec as mvs
from molviewspec.nodes import ComponentExpression
import pandas as pd


# will prob to refactor this later,
# cif parameter only used if color=chains
# so prob would want fix that later
def visualize_macromolecule(structure,
                            cif,

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



def visualize_restraint(structure,

                        start_asym_id: str | int, 
                        start_seq_id: int, 

                        end_asym_id: str | int,
                        end_seq_id: int, 

                        distance: float,
                        restraint_type: str, 

                        representation: str="ball_and_stick",
                        start_atom_id: str="CA", end_atom_id: str="CA",

                        color: str="red", radius: float=0.1, dash_length: float=0.1,
                        label_template="Solved Distance: {{distance}}" + ", Restraint Distance: {restraint_symbol}{distance}",
                        label_color=None,
                        residue_color=None,

                        focus: bool=True,
                        ):

    if start_atom_id == ".":
        start_atom_id = "CA"

    if end_atom_id == ".":
        end_atom_id = "CA"

    hex_color = color_utils.get_color(color)

    if label_color is None:
        label_hex_color = hex_color
    else:
        label_hex_color = color_utils.get_color(label_hex_color)

    if residue_color is None:
        residue_hex_color = hex_color
    else:
        residue_hex_color = color_utils.get_color(residue_hex_color)
        

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

    structure.component(selector=start_residue).representation(type=representation).color(color=residue_hex_color)
    structure.component(selector=end_residue).representation(type=representation).color(color=residue_hex_color)

    res = structure.primitives().distance(
            start=start_atom,
            end=end_atom,
            color=hex_color,
            radius=radius,
            dash_length=dash_length,
            label_template=label_template.format(restraint_symbol=restraint_type_to_symbol(restraint_type), distance=distance),
            label_color=label_hex_color)

    if focus:
        res.focus()






