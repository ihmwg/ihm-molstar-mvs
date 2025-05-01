import pytest
import json
from pathlib import Path
import yaml

import ihm_vis
from ihm_vis import style

##############################################################
# Helper function for 
# testing dict equality
#######################

def dict_equality(dict_A, dict_B):

    dict_A_str = json.dumps(dict_A, sort_keys=True)
    dict_B_str = json.dumps(dict_B, sort_keys=True)

    return dict_A_str == dict_B_str

##############################################################
# test files
############

# Assume tests are always run 
# at top level directory
TEST_ROOT = (Path(".") / "tests").resolve()

TEST_JSON_FILE = TEST_ROOT / "test_data" / "test.json"
assert TEST_JSON_FILE.exists()
TEST_JSON_FILE = str(TEST_JSON_FILE)

TEST_YAML_FILE = TEST_ROOT / "test_data" / "test.yaml"
assert TEST_YAML_FILE.exists()
TEST_YAML_FILE = str(TEST_YAML_FILE)


##############################################################
# user style example
####################

TEST_USER_STYLE = { 

        "macromolecule": {
            "default": {
                 "representation_params": {
                                                "type": "ball_and_stick",
                                           },
                 "opacity_params": {
                                       "opacity": 0.1,
                                   },
                 },
        },

        "component": {
            "default": {
                "opacity_params": {
                                       "opacity": 0.1,
                                  },
                "color_params": {
                                        "color": "blue",
                                },
                },

            "violated": {
                "color_params": {
                                    "color": "red",
                                },
                        },
            },

        "distance": {

                    "default": {
                                 "distance_params": {
                                                      "radius": 0.5,
                                                    },

                               },
                    "violated": {
                                "distance_params": {
                                                     "color": "red",
                                                },

                    },
       },

}


##############################################################
# General Tests
###############

def test_merge():

    base = {
				"top_level_base_key": {

					"second_level_base_key1": {

						"third_level_base_key1": "base_value",
						"third_level_base_key2": "base_value",
						"third_level_base_key3": {
							
									"fourth_level_updater_key": "base_value"
						},

						
						"third_level_base_key4": "base_value",
					},

				"second_level_base_key2": "base_value"
			}
	}


    updater = {
				"top_level_updater_key": "updater_value",
				"top_level_base_key" : {
				   
					"second_level_base_key1": {
					   
						"third_level_base_key1": style.DEFAULT,
						"third_level_base_key2": "updater_value",
						"third_level_base_key3": {
							
							"fourth_level_updater_key": "updater_value"
						},
				   },
				
					"second_level_base_key2": None,
					"second_level_updater_key2": None,
				   }
	}



    expected = {
				"top_level_updater_key": "updater_value",
				"top_level_base_key": {

					"second_level_base_key1": {

						"third_level_base_key1": "base_value",
						"third_level_base_key2": "updater_value",
						"third_level_base_key3": {
							
									"fourth_level_updater_key": "updater_value",
						},

						
						"third_level_base_key4": "base_value",
					},

                    "second_level_base_key2": None,
                    "second_level_updater_key2": None,
			}
	}


    merged = style._merge(base, updater)
    assert dict_equality(merged, expected)




def test_set_style_and_reset_style():

    # Manually clear
    style.USER_STYLE = {}

    style.set_style(TEST_USER_STYLE)
    assert dict_equality(TEST_USER_STYLE, style.USER_STYLE)

    # test reset
    style.reset_style()
    assert len(style.USER_STYLE) == 0



def test_set_style_from_json():

    # Manually clear
    style.USER_STYLE = {}

    style.set_style_from_json(TEST_JSON_FILE)

    with open(TEST_JSON_FILE, "r") as f:
        test_style = json.load(f)

    assert dict_equality(test_style, style.USER_STYLE)

    # test reset
    style.reset_style()
    assert len(style.USER_STYLE) == 0

def test_set_style_from_yaml():

    # Manually clear
    style.USER_STYLE = {}

    style.set_style_from_yaml(TEST_YAML_FILE)

    with open(TEST_YAML_FILE, "r") as f:
        test_style = yaml.safe_load(f)

    assert dict_equality(test_style, style.USER_STYLE)

    # test reset
    style.reset_style()
    assert len(style.USER_STYLE) == 0


##############################################################
# Macromolecule
###############

def test_resolve_style_macromolecule():

    # Set user style
    style.reset_style()
    style.set_style(TEST_USER_STYLE)

    # params
    params = {"representation_params": {
                                       "size_factor": 0.7
                                       },

            "opacity_params": style.DEFAULT,
                                        
                
            "color_params": {
                                    "custom": None,
                            }
             }


    resolved_style = style.resolve_style("macromolecule", params)

    # user set styles
    assert resolved_style["representation_params"]["type"] == "ball_and_stick"
    assert dict_equality(resolved_style["opacity_params"], TEST_USER_STYLE["macromolecule"]["default"]["opacity_params"])

    # param set styles
    assert resolved_style["representation_params"]["size_factor"] == 0.7
    assert resolved_style["color_params"]["custom"] is None

    # default inherited values
    assert resolved_style["color_params"]["color"] == style.DEFAULT_STYLE["macromolecule"]["default"]["color_params"]["color"]



##############################################################
# Components
############

def test_resolve_style_component():

    # Set user style
    style.reset_style()
    style.set_style(TEST_USER_STYLE)

    # params
    params = {"representation_params": {
                                       "size_factor": 0.1
                                       },

            "opacity_params": style.DEFAULT,
                                        
                
            "color_params": {
                                    "custom": None,
                            }
             }


    resolved_style = style.resolve_style("component", params)

    # user set styles
    assert dict_equality(resolved_style["opacity_params"], TEST_USER_STYLE["component"]["default"]["opacity_params"])

    # param set styles
    assert resolved_style["representation_params"]["size_factor"] == 0.1
    assert resolved_style["color_params"]["custom"] is None

    # default inherited values
    assert resolved_style["representation_params"]["type"] == style.DEFAULT_STYLE["component"]["default"]["representation_params"]["type"]


 
def test_resolve_style_component_sub_style():

    # Set user style
    style.reset_style()
    style.set_style(TEST_USER_STYLE)

    # params
    params = {"representation_params": {
                                       "size_factor": 0.1
                                       },

            "opacity_params": style.DEFAULT,
                                        
                
            "color_params": {
                                    "custom": None,
                            }
             }


    resolved_style = style.resolve_style("component", params, "violated")

    # user set styles
    assert dict_equality(resolved_style["opacity_params"], TEST_USER_STYLE["component"]["default"]["opacity_params"])

    # Sub_style override
    assert resolved_style["color_params"]["color"] == "red"

    # param set styles
    assert resolved_style["representation_params"]["size_factor"] == 0.1
    assert resolved_style["color_params"]["custom"] is None

    # default inherited values
    assert resolved_style["representation_params"]["type"] == style.DEFAULT_STYLE["component"]["default"]["representation_params"]["type"]


##############################################################
# Distances
###########

def test_resolve_style_distance():

    # Set user style
    style.reset_style()
    style.set_style(TEST_USER_STYLE)

    # params
    params = {"distance_params": {
                                        "label_size": 1.8,
                                        "label_template": style.DEFAULT
                                },
             }


    resolved_style = style.resolve_style("distance", params)

    # user set styles
    resolved_style["distance_params"]["radius"] == TEST_USER_STYLE["distance"]["default"]["distance_params"]["radius"]

    # param set styles
    assert resolved_style["distance_params"]["label_size"] == 1.8

    # default inherited values
    assert resolved_style["distance_params"]["dash_length"] == style.DEFAULT_STYLE["distance"]["default"]["distance_params"]["dash_length"]
    assert resolved_style["distance_params"]["label_template"] == style.DEFAULT_STYLE["distance"]["default"]["distance_params"]["label_template"]

 
def test_resolve_style_distance_sub_style():

    # Set user style
    style.reset_style()
    style.set_style(TEST_USER_STYLE)

    # params
    params = {"distance_params": {
                                        "label_size": 1.8,
                                        "label_template": style.DEFAULT,
                                },
             }

    resolved_style = style.resolve_style("distance", params, "violated")

    # user set styles
    resolved_style["distance_params"]["radius"] == TEST_USER_STYLE["distance"]["default"]["distance_params"]["radius"]

    # sub_style overrides
    assert resolved_style["distance_params"]["color"] == TEST_USER_STYLE["distance"]["violated"]["distance_params"]["color"]

    # param set styles
    assert resolved_style["distance_params"]["label_size"] == 1.8

    # default inherited values
    assert resolved_style["distance_params"]["dash_length"] == style.DEFAULT_STYLE["distance"]["default"]["distance_params"]["dash_length"]
    assert resolved_style["distance_params"]["label_template"] == style.DEFAULT_STYLE["distance"]["default"]["distance_params"]["label_template"]

