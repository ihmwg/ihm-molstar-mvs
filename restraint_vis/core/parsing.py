"""
Functions for parsing cif files and extracting restraint information
"""

from mmcif.io.PdbxReader import PdbxReader
import pandas as pd
import io
import requests
from urllib.parse import urlparse
from pathlib import Path

def parse_cif_from_file(filepath):
    """
    Parse an mmCIF (macromolecular Crystallographic Information File) file and return its contents.

    Parameters 
    ----------
    filepath : str
        The path to the mmCIF file that will be parsed.

    Returns
    -------
    containers : list
        A list containg the parsed data from the mmCIF file. 
    """
    # list to hold all data
    containers = []
    # open the mmcif file
    with open(filepath, 'r') as ifh:
        # load data into containers list
        PdbxReader(ifh).read(containers)
    return containers

def parse_cif_from_url(url):
    """
    Parse an mmCIF (macromolecular Crystallographic Information File) file from a url and return its contents.

    Parameters 
    ----------
    url : str
        The url of the mmCIF file that will be parsed.

    Returns
    -------
    containers : list
        A list containg the parsed data from the mmCIF file. 
    """
    # list to hold all data
    containers = []
    # get file from URL
    response = requests.get(url)
    # Check if response was successful
    if response.status_code == 200:
        file_content = io.StringIO(response.text)
        # load data into containers list
        PdbxReader(file_content).read(containers)
    else:
        raise ValueError(f"Failed to retrieve file from {url}. Input must be a valid url")
    return containers

def get_atom_coordinates(container, atom_id, comp_id, entity_id, asym_id, seq_id):
    """
    Retrieve the Cartesian coordinates (x, y, z) of an atom based on its identifiers.

    Parameters
    ----------
    container : object
        The container object which contains atom site data.

    atom_id : str
        The identifier for the atom.
    
    comp_id : str
        The component identifer to match for the atom.
    
    entity_id : str
        The entity identifier to match for the atom.
    
    asym_id : str
        The asymmetry identifier to match for the atom.
    
    seq_id : str
        The sequence identifier to match for the atom.
    
    Returns
    -------
    tuple or None
        A tuple (x, y, z) representing the Cartesian coordinates of the atom if a 
        matching atom is found. If no matching atom is found, None is returned.
    """
    # Atom coordinates are in the 'atom_site' category
    atom_data = container.getObj('atom_site')
    # Coordinates from primitive in 'ihm_sphere_obj_site' category
    sphere_data = container.getObj('ihm_sphere_obj_site')
    if atom_data is None and sphere_data is None:
        raise ValueError("Neither 'atom_site' nor 'ihm_sphere_obj_site' categories are available in the CIF file.")
    if atom_data is not None:
        for i in range(atom_data.getRowCount()):
            atom_ids = atom_data.getValue('label_atom_id', i)
            comp_ids = atom_data.getValue('label_comp_id', i)
            entity_ids = atom_data.getValue('label_entity_id', i)
            asym_ids = atom_data.getValue('label_asym_id', i)
            seq_ids = atom_data.getValue('label_seq_id', i)
            coords_x = atom_data.getValue('Cartn_x', i)
            coords_y = atom_data.getValue('Cartn_y', i)
            coords_z = atom_data.getValue('Cartn_z', i)
            # cross reference to ids from cross link
            if (atom_ids == atom_id and
                comp_ids == comp_id and
                entity_ids == entity_id and
                asym_ids == asym_id and
                seq_ids == seq_id):
                x, y, z, = coords_x, coords_y, coords_z
                return (float(x), float(y), float(z)) # match found 
        return None # no match
    elif sphere_data is not None:
        for i in range(sphere_data.getRowCount()):
            entity_ids = sphere_data.getValue('entity_id', i)
            asym_ids = sphere_data.getValue('asym_id', i)
            seq_begin_ids = sphere_data.getValue('seq_id_begin', i)
            seq_end_ids = sphere_data.getValue('seq_id_end', i)
            coords_x = sphere_data.getValue('Cartn_x', i)
            coords_y = sphere_data.getValue('Cartn_y', i)
            coords_z = sphere_data.getValue('Cartn_z', i)
            # cross reference to ids from cross link
            if (entity_ids == entity_id and
                asym_ids == asym_id and
                seq_begin_ids == seq_id or
                seq_end_ids == seq_id):
                x, y, z, = coords_x, coords_y, coords_z
                return (float(x), float(y), float(z)) # match found
        return None # no match

def get_restraints(containers):
    """
    Extract cross-link restraint data from a list of containers and return it as a pandas DataFrame.

    Parameters
    ----------
    containers : list
        A list of containers, each containing objects with cross-link restraint data.

    Returns
    -------
    df : pd.DataFrame
        A pandas DataFrame where each row represents a cross-link restraint, ith columns including:
        'entity_id_1', 'asym_id_1', 'seq_id_1', 'comp_id_1', 'atom_id_1', 
        'entity_id_2', 'asym_id_2', 'seq_id_2', 'comp_id_2', 'atom_id_2', 
        'model_granularity', 'distance_threshold', and 'restraint_type'.
        Rows will be empty if no cross-link restraints found.
    """
    for container in containers:
        restraints = []
        # Extract cross-link restraint data from the container
        cross_link_data = container.getObj('ihm_cross_link_restraint')
        # Fallback: no restraints will return empty dataframe
        if (cross_link_data == None):
            return pd.DataFrame(columns=["entity_id_1", "asym_id_1", "seq_id_1", "comp_id_1", "atom_id_1",
                                         "entity_id_2", "asym_id_2", "seq_id_2", "comp_id_2", "atom_id_2",
                                         "model_granularity", "distance_threshold", "restraint_type",
                                         "atom_id_1_coords", "atom_id_2_coords"])
        # Iterate through each row of cross-link data
        for i in range(cross_link_data.getRowCount()):
            data = {
            "entity_id_1":  cross_link_data.getValue("entity_id_1", i),
            "asym_id_1": cross_link_data.getValue("asym_id_1", i),
            "seq_id_1": cross_link_data.getValue("seq_id_1", i),
            "comp_id_1": cross_link_data.getValue("comp_id_1", i),
            "atom_id_1" : cross_link_data.getValue("atom_id_1", i),

            "entity_id_2": cross_link_data.getValue("entity_id_2", i),
            "asym_id_2": cross_link_data.getValue("asym_id_2", i),
            "seq_id_2": cross_link_data.getValue("seq_id_2", i),
            "comp_id_2": cross_link_data.getValue("comp_id_2", i),
            "atom_id_2" : cross_link_data.getValue("atom_id_2", i),

            "model_granularity": cross_link_data.getValue("model_granularity", i),
            "distance_threshold": float(cross_link_data.getValue("distance_threshold", i)),
            "restraint_type": cross_link_data.getValue("restraint_type", i)
            }
            restraints.append(data)
    # Convert the list of restraints to a pandas DataFrame.
    df = pd.DataFrame(restraints)
    # If atom ids are not specified default to carbon alpha (CA)
    df['atom_id_1'] = df['atom_id_1'].str.replace('.', 'CA', regex=False)
    df['atom_id_2'] = df['atom_id_2'].str.replace('.', 'CA', regex=False)
    # Lists to store retrieved atom coordinates
    atom_id_1_coords = []
    atom_id_2_coords = []
    # Iterate through cross link in dataframe
    for index, row in df.iterrows():
        # Get coordinates for each atom in crosslink
        atom_id_1_coords.append(get_atom_coordinates(container, row['atom_id_1'], row['comp_id_1'], row['entity_id_1'], row['asym_id_1'], row['seq_id_1']))
        atom_id_2_coords.append(get_atom_coordinates(container, row['atom_id_2'], row['comp_id_2'], row['entity_id_2'], row['asym_id_2'], row['seq_id_2']))
    # Add coordinates to dataframe
    df['atom_id_1_coords'] = atom_id_1_coords
    df['atom_id_2_coords'] = atom_id_2_coords
    return df

def get_restraints_dataframe(source):
    """
    Top level wrapper function to take in an mmCIF (from url or file path) and output its restraint DataFrame.

    The function handles input source resolution (remote or local), parses the mmCIF file accordingly, 
    and calls the `get_restraints` method to extract the cross link restraints.

    Parameters 
    ----------
    source : str
        The url or local file path of the mmCIF file that will be parsed.

    Returns
    -------
    restraint_df : pd.DataFrame
        A pandas Aataframe containing crosslink restraint data with defined columns, 
        or an empty data DataFrame with the same structure if no crosslink restraints are found.
    """
    # Psrse the input source to determine if it's a URL or local file path
    parsed = urlparse(source)
    # If source is URL, will call `parse_cif_from_url`
    if parsed.scheme in ('http', 'https'):
        containers = parse_cif_from_url(source)
    # If source is a valid local file path, will call `parse_cif_from_file`
    elif Path(source).exists():
        containers = parse_cif_from_file(source)
    # If the source is neither a valid url of file path, raise an error
    else:
        raise ValueError("Invalid source: must be a valid URL or existing file path.")
    # Extract crosslink restraints from parsed mmCIF containers
    restraint_df = get_restraints(containers)
    return restraint_df