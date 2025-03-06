"""
Misc. utilities
"""

#from typing import
import operator


RESTRAINT_TYPE_TO_SYMBOL = {

        "upper bound": "<",
        "lower bound": ">",

        }

RESTRAINT_TYPE_TO_OPERATOR = {

        "upper bound": operator.lt,
        "lower bound": operator.gt,

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



def restraint_type_to_operator(restraint_type: str):
    """
    restraint type to operator

    raises error if symbol is not known

    Parameters
    ----------
    restraint_type : str, restraint type from mmcif

    Returns
    -------
    symbol: Callable
        corresponding operator
    """

    
    op = RESTRAINT_TYPE_TO_OPERATOR.get(restraint_type, None)
    if op is None:
        raise ValueError("dont know which operator to use for restraint type: " + restraint_type)

    return op




