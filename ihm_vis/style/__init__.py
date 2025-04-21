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
        if isinstance(u, Default):
            continue

        # Not an iterable, uhh
        if not hasattr(u, "items"):
            continue

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


def resolve_style(type, params, sub_style=None):

    defaults = DEFAULT_STYLE.get(type, {}).get("default", {})
    user = USER_STYLE.get(type, {}).get("default", {})

    if sub_style is None:
        return _merge(defaults, user, params)

    defaults_sub_style = DEFAULT_STYLE.get(type, {}).get(sub_style, {})
    user_sub_style = USER_STYLE.get(type, {}).get(sub_style, {})
    return _merge(defaults, defaults_sub_style, user, user_sub_style, params)


class MacromoleculeStyle:
    def __init__(self, selector, representation_params=DEFAULT, color_params=DEFAULT, opacity_params=DEFAULT):

        self.selector = selector

        params = {"representation_params": representation_params,
                  "color_params": color_params,
                  "opacity_params": opacity_params}

        params = resolve_style("macromolecule", params)

        self.representation_params = params["representation_params"]
        self.color_params = params["color_params"]
        self.opacity_params = params["opacity_params"]


    def update(self, representation_params={}, color_params={}, opacity_params={}):
        self.representation_params = _merge(self.representation_params, representation_params)
        self.color_params = _merge(self.color_params, color_params)
        self.opacity_params = _merge(self.opacity_params, opacity_params)


class ComponentStyle:
    def __init__(self, selector, representation_params=DEFAULT, color_params=DEFAULT, opacity_params=DEFAULT, macromolecule_opacity_params=DEFAULT, sub_style="default"):

        self.selector = selector
        self.sub_style = sub_style

        params = {"representation_params": representation_params,
                  "color_params": color_params,
                  "opacity_params": opacity_params,
                  "macromolecule_opacity_params": macromolecule_opacity_params}

        params = resolve_style("component", params, sub_style=self.sub_style)

        self.representation_params = params["representation_params"]
        self.color_params = params["color_params"]
        self.opacity_params = params["opacity_params"]
        self.macromolecule_opacity_params = params["macromolecule_opacity_params"]


    def update(self, representation_params={}, color_params={}, opacity_params={}, macromolecule_opacity_params={}, sub_style=None):

        # if sub_style has changed, need complete reinit
        if sub_style is not None and sub_style != self.sub_style:
            self.__init__(self.selector, representation_params, color_params, opacity_params, macromolecule_opacity_params, sub_style)

        # otherwise just update
        else:
            self.representation_params = _merge(self.representation_params, representation_params)
            self.color_params = _merge(self.color_params, color_params)
            self.opacity_params = _merge(self.opacity_params, opacity_params)
            self.macromolecule_opacity_params = _merge(self.macromolecule_opacity_params, macromolecule_opacity_params)


class DistanceStyle:
    def __init__(self, start_selector, end_selector, distance_params, sub_style="default", **label_keys):

        self.start_selector = start_selector
        self.end_selector = end_selector
        self.sub_style = sub_style
        self.label_keys = label_keys

        params = {"distance_params": distance_params}
        params = resolve_style("distance", params, sub_style=self.sub_style)

        self.distance_params = params["distance_params"]

    def update(self, distance_params={}, sub_style=None, **label_keys):

        # if sub_style has changed, need complete reinit
        if sub_style is not None and sub_style != self.sub_style:
            print("current sub style", self.sub_style, "new sub style", sub_style)
            self.__init__(self.start_selector, self.end_selector, distance_params, sub_style, **label_keys)

        # otherwise just update
        else:
            self.distance_params = _merge(self.distance_params, distance_params)
            self.label_keys = _merge(self.label_keys, label_keys)


class StyleDict:
    def __init__(self):

        self.macromolecule = None
        self.components = {}
        self.distances = {}


    def set_macromolecule_style(self, selector,
                                representation_params=DEFAULT,
                                color_params=DEFAULT,
                                opacity_params=DEFAULT):

        if self.macromolecule is None:
            self.macromolecule = MacromoleculeStyle(selector, 
                                                       representation_params,
                                                       color_params,
                                                       opacity_params)

        else:
            self.macromolecule.update(
                                                    representation_params,
                                                    color_params,
                                                    opacity_params)

        
    def set_component_style(self, selector, 
                              representation_params=DEFAULT, 
                              color_params=DEFAULT,
                              opacity_params=DEFAULT,
                              macromolecule_opacity_params=DEFAULT,
                              sub_style="default"):

        if isinstance(selector, str):
            key = selector

        elif hasattr(selector, "json"):
            key = selector.json()

        else:
            raise ValueError("Unknown selector type")


        if not key in self.components:
            self.components[key] = ComponentStyle(selector, 
                                                       representation_params,
                                                       color_params,
                                                       opacity_params,
                                                       macromolecule_opacity_params,
                                                       sub_style)

        else:
            self.components[key].update(
                                                    representation_params,
                                                    color_params,
                                                    opacity_params,
                                                    macromolecule_opacity_params,
                                                    sub_style)

        # Little weird, but each component has the ability
        # to also affect the macromolecule opacity
        # so here, if a component has set a particular macromolecule opacity
        # make that change refected on the macromolecule
        if self.components[key].macromolecule_opacity_params:
            self.macromolecule.update(opacity_params=self.components[key].macromolecule_opacity_params)


    def set_distance_style(self, start_selector, end_selector, distance_params=DEFAULT, sub_style="default", **label_keys):

        if isinstance(start_selector, str):
            start_key = start_selector

        elif hasattr(start_selector, "json"):
            start_key = start_selector.json()

        else:
            raise ValueError("Unknown selector type")


        if isinstance(end_selector, str):
            end_key = end_selector

        elif hasattr(end_selector, "json"):
            end_key = end_selector.json()

        else:
            raise ValueError("Unknown selector type")

        key = (start_key, end_key)

        if not key in self.distances:
            self.distances[key] = DistanceStyle(start_selector, end_selector,
                                                     distance_params,
                                                     sub_style,
                                                     **label_keys, 
                                                     )

        else:
            self.distances[key].update(
                                          distance_params,
                                          sub_style,
                                          **label_keys,
                                                    )


#######################################

def set_style(kwargs):
    """
    Update default style with user preferences

    The style information is stored as a nested dictionary of {"function_name": {"arg": "default"}}. Therefore, to change the default restraint color to #f0624d you would call set_style({"visualize_restraint": {"color": "#f0624d"}}). Use the get_style function to return the current style dictionary and see availible function and argument names.

    Parameters
    ----------
    kwargs: mappings of "function_name" : {"arg": "user_specified_default"}

    """
    USER_STYLE.update(kwargs)

def set_style_from_json(file: str|Path):
    """
    Update default style with user preferences from json file

    The style information is stored as a nested dictionary of {"function_name": {"arg": "default"}}. Therefore, to change the default restraint color to #f0624d you would call set_style({"visualize_restraint": {"color": "#f0624d"}}). Use the get_style function to return the current style dictionary and see availible function and argument names.

    Parameters
    ----------
    file : json file containing mappings of "function_name" : {"arg": "user_specified_default"}

    """

    with open(file, "r") as f:
        USER_STYLE.update(json.load(f))


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
        USER_STYLE.update(yaml.safe_load(f))


def reset_style():
    """
    Remove all user-specified style preferences and return to defaults defined by restraint_vis/config/defaults.json
    """
    USER_STYLE.clear()

def get_style() -> Dict:
    """
    Return dictionary of current style

    This should not be used to modify the style, only to view. To update, see :func:`restraint_vis.config.set_style`

    Returns
    ______
       style: Dict
           current style dictionary
    """
    return _merge(DEFAULT_STYLE, USER_STYLE) 


