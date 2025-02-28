"""
Misc. utilities
"""

#from typing import


RESTRAINT_TYPE_TO_SYMBOL = {

        "upper bound": "<",
        "lower bound": ">",

        }

def restraint_type_to_symbol(restraint_type: str) -> str:
    """
    restraint type to symbol

    returns '?' if symbol is not known

    Parameters
    ----------
    restraint_type : str, restraint type from mmcif

    Returns
    -------
    symbol: str
        symbol rep
    """

    return RESTRAINT_TYPE_TO_SYMBOL.get(restraint_type, "?")





