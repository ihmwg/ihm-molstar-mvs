import pytest

import pandas as pd
import numpy as np
from ihm_vis.utils import restraint_filters

##############################################################
# Dataframe for testing
#######################

def make_test_df():
    return pd.DataFrame(

        [

            # Different chains, distance 0, compliant
            ("1", "A", "677", "LYS", "CA",   "1", "B", "47", "LYS", "CA",   "by-residue", 35.0, "upper bound", (0, 0, 0), (0, 0, 0)),

            # Same chain, distance 0, compliant
            ("1", "A", "677", "LYS", "CA",   "1", "A", "47", "LYS", "CA",   "by-residue", 35.0, "upper bound", (0, 0, 0), (0, 0, 0)),

            # Same chain, distance 100, violated
            ("1", "A", "677", "LYS", "CA",   "1", "A", "47", "LYS", "CA",   "by-residue", 35.0, "upper bound", (100, 0, 0), (0, 0, 0)),

            # Same chain, distance 0, both close to 100, 100, 100, compliant
            ("1", "A", "677", "LYS", "CA",   "1", "A", "47", "LYS", "CA",   "by-residue", 35.0, "upper bound", (100, 100, 100), (100, 100, 100)),

            # Same chain, distance far, one close to 100, 100, 100, violated
            ("1", "A", "677", "LYS", "CA",   "1", "A", "47", "LYS", "CA",   "by-residue", 35.0, "upper bound", (100, 100, 100), (-100, -100, -100)),


        ], columns = ["entity_id_1",
                      "asym_id_1",
                      "seq_id_1",
                      "comp_id_1",
                      "atom_id_1",

                      "entity_id_2",
                      "asym_id_2",
                      "seq_id_2",
                      "comp_id_2",
                      "atom_id_2",

                      "model_granularity",
                      "distance_threshold",
                      "restraint_type",

                      "atom_id_1_coords",
                      "atom_id_2_coords",
                      ]
    )


##############################################################
# Get solved distance
#####################

def test_get_solved_distance():

    df = make_test_df()
    df = restraint_filters.get_solved_distance(df)

    assert np.allclose(df["solved_distance"].values, [0, 0, 100, 0, 346.4101615])
    assert np.all(df["compliant"].values == [True, True, False, True, False])



##############################################################
# Each filter func
##################

def test_first_n():

    df_A = make_test_df()
    df_A = restraint_filters.BUILTIN_FILTER_FUNCS["first_n"](df_A, n=1)

    assert len(df_A) == 1

    df_B = make_test_df()
    df_B = restraint_filters.BUILTIN_FILTER_FUNCS["first_n"](df_B, n=3)

    print(df_B)

    assert len(df_B) == 3


def test_random_sample():

    # not sure if I should test for randomness
    # seems overkill, I think really this is testing
    # whether the builtin is available bc otherwise 
    # we are just testing pandas itself...

    df_A = make_test_df()
    df_A = restraint_filters.BUILTIN_FILTER_FUNCS["random_sample"](df_A, n=1)

    assert len(df_A) == 1

    df_B = make_test_df()
    df_B = restraint_filters.BUILTIN_FILTER_FUNCS["random_sample"](df_B, n=3)

    assert len(df_B) == 3

def test_across_chains():

    df = make_test_df()
    df = restraint_filters.BUILTIN_FILTER_FUNCS["across_chains"](df)

    # should only retain the second row
    assert len(df.index) == 1
    assert df.index[0] == 0


def test_violated():

    df = make_test_df()
    df = restraint_filters.BUILTIN_FILTER_FUNCS["violated"](df)

    assert len(df) == 2
    assert np.all(df.index.values == [2, 4])


def test_compliant():

    df = make_test_df()
    df = restraint_filters.BUILTIN_FILTER_FUNCS["compliant"](df)

    assert len(df) == 3
    assert np.all(df.index.values == [0, 1, 3])


def test_within_X_of():

    df = make_test_df()
    df = restraint_filters.BUILTIN_FILTER_FUNCS["within_X_of"](df, X=0.1, position=(100, 100, 100))

    assert len(df) == 2
    assert np.all(df.index.values == [3, 4])






