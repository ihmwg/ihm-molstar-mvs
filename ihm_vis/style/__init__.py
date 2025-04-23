"""
Configuration management utilities for visualization styles.

This module provides a system for defining, updating, and resolving style configurations 
for macromolecular and component visualizations, including distances between components. 
It supports hierarchical merging of default and user-specified styles from dictionaries, 
JSON files, and YAML files. The system allows users to apply fine-grained control over 
representation, color, and opacity parameters.

Configuration Description
--------
The style system is built around “types” of visuals — macromolecule, component, and distance — each of which can have multiple sub_styles. The system is hierarchical and stored in a series of nested dictionaries of the form:

.. code-block:: python

    {
    "<visual_type>": {
        "<sub_style_name>": {   # All 
            # parameter groups
            "representation_params": { … },
            "color_params":          { … },
            "opacity_params":        { … },
        },
        # (optionally) other sub_styles for this visual_type
    },
    # next visual_type…
    }



The top level keys are for the "macromolecule", the main structure you are visualizing, the "component"s, which are the individual residues within a restraint, and "distance"s, the label primities drawn between the individual "component"s of the restraint.

Why "component" instead of "residue"? We wanted to be more generic in case future IHM data will have restraints of different types. Though, the current functionality is geared toward cross-linking experiments and generally when we refer to "component"s in this documentation you can internalize that as "residue".

The second level of keys defines a "sub_style" within that visual type category. There is always a base sub_style key named "default". This sets the default visualization of that type - for example the following sets the default visualization of any restraint residues to a blue ball and stick.

.. code-block:: python

    "component": {
        "default": {
            "representation_params": {"type": "ball_and_stick"}
            "color_params": {"color": "blue"}
        }
    }


This means any component we visualize will inherit this ball and stick representation unless explicitly overridden. However, sometimes we may want to have some components styled differently than others. To achieve this, the user can either pass specific parameters to the `IHM_Builder.set_single_restraint_style`, or can specify additional "sub_style" as follows:


.. code-block:: python

    "component": {
        "default": {
            "representation_params": {"type": "ball_and_stick"},
            "color_params": {"color": "blue"}
        },

        "violated": {
            "color_params": {"color": "red"}
        }

        "compliant": {
            "color_params": {"color": "green"}
        }
    }


Here, we have defined two sub_styles called "violated" and "compliant". The user can now pass `sub_style="violated"` 




Main elements for developers
----------------------------
- DEFAULT_STYLE, USER_STYLE are global style dictionaries
- _merge function is used to 

Intended for use with restraint-based macromolecular visualization pipelines, where consistent, reusable,
and customizable styling is essential.

Example
"""



from pathlib import Path
from collections.abc import Mapping
from typing import Dict, Any, Optional
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
    """
    Recursively merge one or more dictionaries into a base dictionary.

    Handles nested mappings and respects special `Default` placeholder objects 
    by omitting or initializing values accordingly.

    Parameters
    ----------
    base : dict
        Base dictionary to merge into.
    *args : dict
        One or more dictionaries to merge into the base.

    Returns
    -------
    dict
        A new dictionary representing the merged result.
    """

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


def resolve_style(type: str, params: Dict[str, Dict[Any, Any]], sub_style: Optional[str]=None):
    """
    Resolve a style dictionary by merging defaults, user settings, and runtime parameters.

    Supports merging in a style hierarchy: global defaults → sub-style defaults → runtime parameters.

    Parameters
    ----------
    type : str
        The style type (e.g., "macromolecule", "component", "distance").
    params : Dict[str, Dict[Any, Any]]
        Runtime parameters to merge into the resolved style.
    sub_style : str, optional
        Optional sub-style to apply if defined in the style configuration.

    Returns
    -------
    dict
        The fully resolved style dictionary.
    """

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
    Update the global user style dictionary with new preferences.

    See :ref:`ihm_vis.style.sub_style_modes` module page for detailed structure of style dict hierachy

    Parameters
    ----------
    kwargs : dict
        Mapping of style types and their associated user-defined parameter overrides.
        For example:
        {"macromolecule": {"default": {"color_params": {"color": "blue"}}}}

    """
    USER_STYLE.update(kwargs)

def set_style_from_json(file: str|Path):
    """
    Load and apply style preferences from a JSON file.

    See :ref:`ihm_vis.style.sub_style_modes` module page for detailed structure of style dict hierachy

    Parameters
    ----------
    file : str or Path
        Path to a JSON file containing user style definitions.
    """

    with open(file, "r") as f:
        USER_STYLE.update(json.load(f))


def set_style_from_yaml(file: str|Path):
    """
    Load and apply style preferences from a YAML file.

    See :ref:`ihm_vis.style.sub_style_modes` module page for detailed structure of style dict hierachy

    Parameters
    ----------
    file : str or Path
        Path to a YAML file containing user style definitions.
    """

    with open(file, "r") as f:
        USER_STYLE.update(yaml.safe_load(f))


def set_style_from_file(file: str|Path):

    file = Path(file)
    if file.suffix == ".yaml" or file.suffix == ".yml":
        set_style_from_yaml(file)

    elif file.suffix == ".json":
        set_style_from_json(file)

    else:
        raise ValueError("Please provide either a .yaml, .yml, or .json file with the requested style")


def reset_style():
    """
    Clear all user-defined styles and revert to the default configuration.

    This resets the `USER_STYLE` dictionary, leaving only the defaults defined by th package itself. 
    """
    USER_STYLE.clear()

def get_style() -> Dict:
    """
    Get the merged style configuration including defaults and user overrides.

    This is a copy and should not be used to modify the configuration directly.

    Returns
    -------
    dict
        The current effective style configuration.
    """
    return _merge(DEFAULT_STYLE, USER_STYLE) 


