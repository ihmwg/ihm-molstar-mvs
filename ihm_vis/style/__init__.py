"""
Helpers for setting global configuration
"""

from pathlib import Path
import inspect
from collections.abc import Mapping
from typing import Dict
import functools
import json
import yaml

# Default placeholder to differentiate from None value
class Default:
    def __init__(self):
        pass

DEFAULT = Default()

USER_STYLE = {}
DEFAULT_STYLE = {}

# Load from defaults file
_config_dir = Path(__file__).parent
_default_config_file = _config_dir / "defaults.json"
with open(_default_config_file, "r") as f:
    DEFAULT_STYLE.update(json.load(f))

def _merge(base: Dict, *args: Dict):
    d = base.copy()

    for u in args:
        for k,v in u.items():
            
            if isinstance(v, Default):
                if not k in d:
                    d[k] = {}

            elif v is None:
                d[k] = None

            elif isinstance(v, Mapping):
                if k in d:
                    d[k] = _merge(d[k], v)
                else:
                    d[k] = v
            else:
                d[k] = v

    # Replace Default
    # with {} for all instances
    # on this level
    return d


def apply_style_defaults(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        # func inputs
        sig = inspect.signature(func)
        bound_args = sig.bind(*args, **kwargs)
        bound_args.apply_defaults()
        bound_args = bound_args.arguments

        # apply func inputs
        # on top of current style

        if "sub_style" in bound_args:
            sub_style = bound_args["sub_style"] 
        else:
            sub_style = "default"
            
        defaults = DEFAULT_STYLE.get(func.__name__, {}).get("default", {})
        defaults_sub_style = DEFAULT_STYLE.get(func.__name__, {}).get(sub_style, {})

        user = USER_STYLE.get(func.__name__, {}).get("default", {})
        user_sub_style = USER_STYLE.get(func.__name__, {}).get(sub_style, {})
            
        bound_args = _merge(defaults, defaults_sub_style, user, user_sub_style, bound_args)

        return func(**bound_args)
    return wrapper


#######################################

def set_style(kwargs):
    """
    Update default style with user preferences

    The style information is stored as a nested dictionary of {"function_name": {"arg": "default"}}. Therefore, to change the default restraint color to #f0624d you would call set_style({"visualize_restraint": {"color": "#f0624d"}}). Use the get_style function to return the current style dictionary and see availible function and argument names.

    Parameters
    ----------
    kwargs: mappings of "function_name" : {"arg": "user_specified_default"}

    """
    _style.update(kwargs)

def set_style_from_json(file: str|Path):
    """
    Update default style with user preferences from json file

    The style information is stored as a nested dictionary of {"function_name": {"arg": "default"}}. Therefore, to change the default restraint color to #f0624d you would call set_style({"visualize_restraint": {"color": "#f0624d"}}). Use the get_style function to return the current style dictionary and see availible function and argument names.

    Parameters
    ----------
    file : json file containing mappings of "function_name" : {"arg": "user_specified_default"}

    """

    with open(file, "r") as f:
        _style.update(json.load(f))


def set_style_from_yaml(file: str|Path):
    """
    Update default style with user preferences from yaml file

    The style information is stored as a nested dictionary of {"function_name": {"arg": "default"}}. Therefore, to change the default restraint color to #f0624d you can provide a yaml file specifying the analogous dictionary structure. For example, th efollowing yaml would set the default color for the visualize_restraint function:

    visualize_restraint:
      color: "#f0624d"
     
     Use the get_style function to return the current style dictionary and see availible function and argument names.

    Parameters
    ----------
    file : yaml file containing mappings of "function_name" : {"arg": "user_specified_default"}

    """

    with open(file, "r") as f:
        _style.update(yaml.safe_load(f))


def reset_style():
    """
    Remove all user-specified style preferences and return to defaults defined by restraint_vis/config/defaults.json
    """
    _style.clear()

def get_style() -> Dict:
    """
    Return dictionary of current style

    This should not be used to modify the style, only to view. To update, see :func:`restraint_vis.config.set_style`

    Returns
    ______
       style: Dict
           current style dictionary
    """
    return _merge(DEFAULT_STYLE, _style) 


