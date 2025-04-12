ihm-molstar-mvs
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/restraint_vis/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/restraint_vis/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ihm-molstar-mvs/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ihm-molstar-mvs/branch/main)


This repository contains utilities that support data-driven visualization of integrative structures through Molstar using MolViewSpec.

# Testing Instructions 4-12-25
We are currently undergoing a refactor of the code (see [refactoring branch](https://github.com/ihmwg/ihm-molstar-mvs/tree/reorganize_codebase)); therefore, the package documentation is being put on pause! 

To demonstrate the current functionality, we have an example notebook (`devtools/jupyter/Quick_start_example.ipynb`) and a cli (`cli/rcsb_script.py`).

### Example Notebook
To run the example notebook, first clone the repo and switch to this branch
```
git clone https://github.com/ihmwg/ihm-molstar-mvs.git
cd ihm-molstar-mvs
git switch main_from_berkeley
```

Create and activate a fresh conda environment 
```
conda create -n ihm_vis python=3.12
conda activate ihm_vis
```

Pip install with dev dependencies (jupyter)
```
pip install .[dev]
```

Start the jupyter server
```
cd devtools/jupyter
jupyter notebook
```

Connect to the jupyter session with the link output from the `jupyter notebook` command and select the `Quick_start_example.ipynb` notebook. 

From there, the notebook contains (minimal) descriptions of each input and functionality within the notebook, which we will expand on in the full docs.

Running the notebook with no modifications will visualize `9a3v.cif` with only restraints that go across chains, styled by whether they are violated or compliant and output the file `9a3v_inter-chain_restraints.mvsj` for viewing in Mol*.

### CLI


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.10.
