"""
Utilities for producing hex color codes
"""

from typing import List

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex, CSS4_COLORS
import re


HEX_COLOR_PAT = re.compile(r"#?([0-9A-Fa-f]{3}|[0-9A-Fa-f]{6})$")

def is_hex_color(maybe_color: str) -> bool:
    """
    Wrapper for colors from matplotlib

    Parameters
    ----------
    maybe_color : str, string that might be a hex color

    Returns
    -------
    is_color : bool
        whether the string was a valid hex color
    """

    return bool(HEX_COLOR_PAT.search(maybe_color))

def get_color(color: str) -> str:
    """
    Wrapper for colors from matplotlib

    Parameters
    ----------
    color_name: str, color name to retrive or hex itself

    Returns
    -------
    color : str
        hex color

    Raises
    ______
        ValueError: color_name not in matplotlibs CSS4 colors
    """

    if is_hex_color(color):
        return color

    color_hex = CSS4_COLORS.get(color, None)

    if color_hex is None:
        raise ValueError(f"The provided color ({color}) wasn't a hex color code, and also wasn't found in the CSS4 colors list. Pleaes either provide a valid hex color code or chooose from:\n" + "\n".join(list(CSS4_COLORS.keys())))

    return color_hex

def get_n_colors(n: int, colormap_name: str="Set3") -> List[str]:
    """
    Produce n hex colors from the requested colormap

    Parameters
    ----------
    n: int, How many unique colors to retrieve 

    Returns
    -------
    colors : List[str]
        n hex colors
    """

    cmap = plt.get_cmap(colormap_name)
    return [to_hex(cmap(i)) for i in np.linspace(0, 1, n)]
