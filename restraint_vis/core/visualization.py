"""
Underlying functions for visualizing with molviewspec
"""

from typing import Dict, Optional

from restraint_vis.utils import color_utils
from restraint_vis.utils.general import RESTRAINT_TYPE_TO_SYMBOL, restraint_type_to_symbol
from restraint_vis.config import apply_style_defaults, DEFAULT

import molviewspec as mvs
from molviewspec.nodes import ComponentExpression
import pandas as pd
import json
from pathlib import Path



# will prob to refactor this later,
# cif parameter only used if color=chains
# so prob would want fix that later
@apply_style_defaults
def visualize_macromolecule(structure, # molviewspec scene builder
                            cif,       # cif file info

                            color_params: Optional[Dict[str, str]]=DEFAULT,
                            representation_params: Optional[Dict[str, str]]=DEFAULT,
                            sub_style: str="default"
                            ): 
    """
    Visualize the macromolecule

    Parameters
    ----------
    structure : mvs.Structure, The structure to visualize
    color:  the color to make the macromolecule

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
    structure.component(selector="polymer").representation(**representation_params).color(**color_params)


# TODO: feels clunky to use, think about improvements
# do we need some abstraction/ shorthand or are we fine
# with a large argument list
@apply_style_defaults
def visualize_restraint(structure,

                        start_asym_id: str|int, 
                        start_seq_id: int, 

                        end_asym_id: str|int,
                        end_seq_id: int, 

                        distance: float,
                        restraint_type: str, 

                        start_atom_id: str="CA", end_atom_id: str="CA",
                        representation_params: Optional[Dict[str, str]]=DEFAULT,
                        color_params: Optional[Dict[str, str]]=DEFAULT,
                        distance_params: Optional[Dict[str, str]]=DEFAULT,
                        tube_params: Optional[Dict[str, str]]=DEFAULT,

                        sub_style="default",

                        focus: Optional[bool]=True,
                        ):

    """
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
    start_component.representation(**representation_params).color(**color_params)

    end_component = structure.component(selector=end_residue)
    end_component.representation(**representation_params).color(**color_params)

    if distance_params is not None:
        if "tooltip" in distance_params:
            distance_params["tooltip"] = distance_params["tooltip"].format(restraint_type_symbol=restraint_type_to_symbol(restraint_type), distance=distance)

        res = structure.primitives().distance(
                start=start_atom,
                end=end_atom,
                **distance_params)

    if tube_params is not None:
        if "tooltip" in tube_params:
            tube_params["tooltip"] = tube_params["tooltip"].format(restraint_type_symbol=restraint_type_to_symbol(restraint_type), distance=distance)
        res = structure.primitives().tube(
                start=start_atom,
                end=end_atom,
                **tube_params)

    if focus:
        res.focus()


