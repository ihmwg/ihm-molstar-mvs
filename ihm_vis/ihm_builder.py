"""
Main interface for package users to visualize IHM restraint data
"""

from typing import Optional, Tuple, Dict, Callable, Any

import molviewspec as mvs
from molviewspec.nodes import ComponentExpression
from molviewspec.mvsx_converter import mvsj_to_mvsx
from pathlib import Path
from mmcif.io import PdbxReader
import pandas as pd
from urllib.parse import urlparse
import requests
import io

from ihm_vis.style import DEFAULT, StyleDict, ComponentStyle, DistanceStyle
from ihm_vis.utils import restraint_type_to_symbol
from ihm_vis.utils.local_file import LocalFile

class IHM_Builder:

    restraint_df_schema = {"entity_id_1": "str",
                           "asym_id_1": "str", 
                           "seq_id_1": "str",
                           "comp_id_1": "str",
                           "atom_id_1": "str",
                           "entity_id_2": "str",
                           "asym_id_2": "str",
                           "seq_id_2": "str",
                           "comp_id_2": "str",
                           "atom_id_2": "str",
                           "model_granularity": "str",
                           "distance_threshold": "float",
                           "restraint_type": "str",
                           }

    def __init__(self, source: str|Path, structure_index: int=0, format: str="mmcif", macromolecule_selector: str="polymer", port: int=8003):
        """
        Initialize an IHM_Builder which loads a protein structure and prepares
        internal state for visualizing restraint data

        Parameters
        ----------
        source : str or Path
            URL or local path to an mmCIF file containing the target structure.
        structure_index : int, optional
            Index of the model/structure to extract from the mmCIF file (default is 0).
        format : str, optional
            Format string passed to MolViewSpec parser (default is "mmcif").
        macromolecule_selector : str, optional
            Selector used to identify macromolecule in within the mmCIF file (default is "polymer").
        port : int, optional
            Local port to serve files if `source` is a local path (default is 8003).

        Raises
        ------
        ValueError
            If `source` is neither a valid HTTP(S) URL nor an existing file path.
        """

        self.source = source
        self.format = format
        self.structure_index = structure_index
        self.macromolecule_selector = macromolecule_selector

        # Depending on source (url vs local file)
        # read in the cif file
        if urlparse(source).scheme in ("http", "https"):
            self.source_type = "url"

            self.local_file = None
            self.url = source
            self.cif = self.read_cif_url(source, self.structure_index)

        elif Path(source).exists():
            self.source_type = "file"

            self.local_file = LocalFile(source, port)
            self.url = self.local_file.url
            self.cif = self.read_cif_file(source, self.structure_index)

        else:
            raise ValueError("Invalid source: must either be a valid URL or a local file path")


        # MVS does not currently support "updating" a representation
        # while there are other tricks to back-track to a previous state,
        # it is simpler to keep an intermediate state here
        # and only call into MVS once we are ready to write the mvsj
        self.state = StyleDict()
        self.state.set_macromolecule_style(selector=self.macromolecule_selector)

        # Setup restraint_df
        cols = {colname : pd.Series(dtype=t) for colname, t in self.restraint_df_schema.items()}
        self.restraint_df = pd.DataFrame(cols)


    ###############################################################################################
    # General Utility functions
    ###########################

    @classmethod
    def read_cif_file(cls, file_path: str|Path, structure_index: int=0) -> PdbxReader.DataContainer:
        """
        Read and parse an mmCIF file from disk into a container object.

        Parameters
        ----------
        file_path : str or Path
            Path to the local mmCIF file.
        structure_index : int, optional
            Index of the desired model in multi‐model mmCIFs (default is 0).

        Returns
        -------
        container : PdbxReader.DataContainer
            Parsed mmCIF container corresponding to the requested model.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """

        containers = []

        with open(file_path, "r") as f:
            PdbxReader.PdbxReader(f).read(containers)

        return containers[structure_index]

    @classmethod
    def read_cif_url(cls, url: str, structure_index: int=0) -> PdbxReader.DataContainer:
        """
        Fetch an mmCIF file from a URL and parse it into a container.

        Parameters
        ----------
        url : str
            HTTP or HTTPS URL pointing to an mmCIF file.
        structure_index : int, optional
            Index of the desired model in the fetched mmCIF (default is 0).

        Returns
        -------
        container : PdbxReader.DataContainer
            Parsed mmCIF container corresponding to the requested model.

        Raises
        ------
        ValueError
            If the HTTP request does not return status code 200.
        """

        containers = []

        response = requests.get(url)
        if response.status_code != 200:
            raise ValueError(f"Failed to retrieve file from {url}. Input must be a valid url")

        f = io.StringIO(response.text)
        PdbxReader(f).read(containers)

        return containers[structure_index]


    @classmethod
    def get_atom_coordinates(cls, container: PdbxReader.DataContainer, atom_id: str, comp_id: str, entity_id: str, asym_id: str, seq_id: str) -> Optional[Tuple[float, float, float]]:
        """
        Retrieve the Cartesian (x, y, z) coordinates of a specific atom.

        Searches both the standard `atom_site` and, if absent, the
        `ihm_sphere_obj_site` categories in the mmCIF container.

        Parameters
        ----------
        container : PdbxReader.DataContainer
            The parsed mmCIF container from which to fetch coordinates.
        atom_id : str
            Label for the atom (e.g., "CA").
        comp_id : str
            Component identifier (residue name) to match.
        entity_id : str
            Entity identifier to match.
        asym_id : str
            Asymmetry (chain) identifier to match.
        seq_id : str
            Sequence number (label_seq_id) to match.

        Returns
        -------
        coords : tuple of float or None
            If found, returns (x, y, z) coordinates as floats; otherwise, None.
        """

        # Atom coordinates are in the 'atom_site' category
        atom_data = container.getObj('atom_site')
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

        atom_data = container.getObj('ihm_sphere_obj_site')
        if atom_data is not None: 
            # Coordinates from primitive
            for i in range(atom_data.getRowCount()):
                entity_ids = atom_data.getValue('entity_id', i)
                asym_ids = atom_data.getValue('asym_id', i)
                seq_begin_ids = atom_data.getValue('seq_id_begin', i)
                seq_end_ids = atom_data.getValue('seq_id_end', i)
                coords_x = atom_data.getValue('Cartn_x', i)
                coords_y = atom_data.getValue('Cartn_y', i)
                coords_z = atom_data.getValue('Cartn_z', i)
                # cross reference to ids from cross link
                if (entity_ids == entity_id and
                    asym_ids == asym_id and
                    seq_begin_ids == seq_id or
                    seq_end_ids == seq_id):
                    x, y, z, = coords_x, coords_y, coords_z
                    return (float(x), float(y), float(z)) # match found
            return None # no match


    ###############################################################################################
    # Parsing different restraint types here #
    ##########################################

    def get_cross_links(self) -> pd.DataFrame:
        """
         Extract cross-link restraints from the loaded mmCIF and store in restraints_df attribute, returns a view.

        Parses the `ihm_cross_link_restraint` category, fills in default
        atom IDs ("CA") where missing, and appends coordinate columns.

        Returns
        -------
        df : pandas.DataFrame
            Columns: ['entity_id_1', 'asym_id_1', 'seq_id_1', 'comp_id_1',
            'atom_id_1', 'entity_id_2', 'asym_id_2', 'seq_id_2',
            'comp_id_2', 'atom_id_2', 'model_granularity',
            'distance_threshold', 'restraint_type', 'atom_id_1_coords',
            'atom_id_2_coords'].
        """

        # Extract cross-link restraint data from the container
        cross_link_data = self.cif.getObj('ihm_cross_link_restraint')
        if cross_link_data is None:
            cols = {colname : pd.Series(dtype=t) for colname, t in self.restraint_df_schema.items()}
            return pd.DataFrame(cols)

        # Iterate through each row of cross-link data
        restraints = []
        for i in range(cross_link_data.getRowCount()):
            restraints.append((
                cross_link_data.getValue("entity_id_1", i),
                cross_link_data.getValue("asym_id_1", i),
                cross_link_data.getValue("seq_id_1", i),
                cross_link_data.getValue("comp_id_1", i),
                cross_link_data.getValue("atom_id_1", i),
                cross_link_data.getValue("entity_id_2", i),
                cross_link_data.getValue("asym_id_2", i),
                cross_link_data.getValue("seq_id_2", i),
                cross_link_data.getValue("comp_id_2", i),
                cross_link_data.getValue("atom_id_2", i),
                cross_link_data.getValue("model_granularity", i),
                float(cross_link_data.getValue("distance_threshold", i)),
                cross_link_data.getValue("restraint_type", i)))


        #the list of restraints to a pandas DataFrame.
        restraint_df = pd.DataFrame(restraints, columns=self.restraint_df_schema.keys())
        # If atom ids are not specified default to carbon alpha (CA)
        restraint_df['atom_id_1'] = restraint_df['atom_id_1'].str.replace('.', 'CA', regex=False)
        restraint_df['atom_id_2'] = restraint_df['atom_id_2'].str.replace('.', 'CA', regex=False)
        # Lists to store retrieved atom coordinates
        atom_id_1_coords = []
        atom_id_2_coords = []
        # Iterate through cross link in dataframe
        for index, row in restraint_df.iterrows():
            # Get coordinates for each atom in crosslink
            atom_id_1_coords.append(self.get_atom_coordinates(self.cif, row['atom_id_1'], row['comp_id_1'], row['entity_id_1'], row['asym_id_1'], row['seq_id_1']))
            atom_id_2_coords.append(self.get_atom_coordinates(self.cif, row['atom_id_2'], row['comp_id_2'], row['entity_id_2'], row['asym_id_2'], row['seq_id_2']))
        # Add coordinates to dataframe
        restraint_df['atom_id_1_coords'] = atom_id_1_coords
        restraint_df['atom_id_2_coords'] = atom_id_2_coords

        
        self.restraint_df = pd.concat((self.restraint_df, restraint_df)).drop_duplicates()
        return restraint_df 

    
    # Master method to parse all restraint types #
    ##############################################
    def get_all_restraints(self) -> pd.DataFrame:
        """
        Parse all supported restraint types into the restrant_df attribute, returns a view.

        Currently only implements cross-link restraints; additional sources
        (e.g., NMR) may be added in future.

        Returns
        -------
        df : pandas.DataFrame
            Combined DataFrame of all parsed restraints.
        """

        # Call each restraint type
        self.get_cross_links()  # cross_linking

        # TODO
        # FUTURE METHODS
        #self.get_nrm_restraints() # nmr

        return self.restraint_df


    ##############################################################################################
    # Filtering restraints
    ######################

    def filter_restraints(self, filter_func: str|Callable, **kwargs) -> pd.DataFrame:
        """
        Apply a filtering function to the internal restraint DataFrame.

        Parameters
        ----------
        filter_func : str or callable
            If str, must be a function name from :ref:`ihm_vis.utils.restraint_filters`
            If callable, should accept (df: DataFrame, **kwargs) and return a filtered DataFrame.
        **kwargs
            Additional arguments to pass to the filter function.

        Returns
        -------
        df : pandas.DataFrame
            The filtered restraint DataFrame.

        Raises
        ------
        ValueError
            If `filter_func` is a string not found in built-in filters.
        """

        # avoid circular imports
        from ihm_vis.utils.restraint_filters import BUILTIN_FILTER_FUNCS

        if isinstance(filter_func, str):
            _filter_func = BUILTIN_FILTER_FUNCS.get(filter_func, None)

 
            if _filter_func is None:
                raise ValueError(f"The requested built-in filter_func ({filter_func}) could not be found. See ihm_vis.filters for availible functions")

        else:
            _filter_func = filter_func

        self.restraint_df = _filter_func(self.restraint_df, **kwargs)
        return self.restraint_df



    ###############################################################################################
    # Visualization Functions
    #########################

    def set_macromolecule_style(self, 
                                representation_params: Optional[dict[str, str]]=DEFAULT, 
                                color_params: Optional[dict[str, str]]=DEFAULT, 
                                opacity_params: Optional[dict[str, str]]=DEFAULT):
        """
        Configure the style of the macromolecule.

        Parameters
        ----------
        representation_params : dict, optional
            Representation settings passed through to MolViewSpec (e.g., {"type": "cartoon"}).
        color_params : dict, optional
            Color settings for the macromolecule (e.g., {"color": "blue"}).
        opacity_params : dict, optional
            Opacity settings for the macromolecule (e.g., {"opacity": 0.5}).
        """

        self.state.set_macromolecule_state(selector=self.macromolecule_selector,
                                representation_params=representation_params,
                                color_params=color_params,
                                opacity_params=opacity_params)

    def set_single_restraint_style(self,

                            start_asym_id: str|int, 
                            start_seq_id: int, 

                            end_asym_id: str|int,
                            end_seq_id: int, 

                            distance: float,
                            restraint_type: str, 

                            start_atom_id: str="CA", end_atom_id: str="CA",

                            representation_params: Optional[Dict[str, str]]=DEFAULT,
                            color_params: Optional[Dict[str, str]]=DEFAULT,
                            opacity_params: Optional[Dict[str, str]]=DEFAULT,
                            distance_params: Optional[Dict[str, str]]=DEFAULT,
                            label_keys={},

                            sub_style="default",

                            focus: Optional[bool]=False,
                            macromolecule_opacity_params: Optional[Dict[str, str]]=DEFAULT):


        """
        Add styling for a single distance restraint between two residues.

        Parameters
        ----------
        start_asym_id : str or int
            Chain identifier of the first residue.
        start_seq_id : int
            Sequence number of the first residue.
        end_asym_id : str or int
            Chain identifier of the second residue.
        end_seq_id : int
            Sequence number of the second residue.
        distance : float
            Experimentally determined distance threshold of this restraint.
        restraint_type : str
            Restraint operator/type (e.g., 'less_than', 'equal').
        start_atom_id : str, optional
            Atom label in the first residue (default 'CA').
        end_atom_id : str, optional
            Atom label in the second residue (default 'CA').
        representation_params : dict, optional
            Visual representation parameters for the two residues (e.g., {"type": "ball_and_stick"}).
        color_params : dict, optional
            Color specification for the residue representations (e.g., {"color": "red"}).
        opacity_params : dict, optional
            Opacity for the residue representations (e.g., {"opacity": 0.5}).
        distance_params : dict, optional
            Parameters for drawing the distance primitive (e.g., {"radius": 0.5}).
        label_keys : dict, optional
            Keys/values for formatting the templated label/tooltip of the distance primitive (e.g., distance value).
        sub_style : str, optional
            Named sub_style to apply (default 'default'). See :ref:`ihm_vis.style.sub_style_modes` for more information.
        focus : bool, optional
            If True, zoom camera to this restraint (default False).
        macromolecule_opacity_params : dict, optional
            Adjust macromolecule opacity to highlight this restraint.
        """

        start_residue = ComponentExpression(label_asym_id=start_asym_id,
                                           beg_label_seq_id=start_seq_id,
                                           end_label_seq_id=start_seq_id)

        end_residue = ComponentExpression(label_asym_id=end_asym_id,
                                           beg_label_seq_id=end_seq_id,
                                           end_label_seq_id=end_seq_id)

        start_atom = ComponentExpression(label_asym_id=start_asym_id,
                                           beg_label_seq_id=start_seq_id,
                                           end_label_seq_id=start_seq_id,
                                           label_atom_id=start_atom_id)

        end_atom = ComponentExpression(label_asym_id=end_asym_id,
                                           beg_label_seq_id=end_seq_id,
                                           end_label_seq_id=end_seq_id,
                                           label_atom_id=end_atom_id)


        self.state.set_component_style(
                                selector=start_residue,
                                representation_params=representation_params,
                                color_params=color_params,
                                opacity_params=opacity_params,
                                sub_style=sub_style)

        self.state.set_component_style(
                                selector=end_residue,
                                representation_params=representation_params,
                                color_params=color_params,
                                opacity_params=opacity_params,
                                sub_style=sub_style)


        # Provide the distance, restraint type, and threshold symbol
        # as automatic label_keys
        if not "distance" in label_keys:
            label_keys["distance"] = distance

        if not "restraint_type" in label_keys:
            label_keys["restraint_type"] = restraint_type

        if not "restraint_type_symbol" in label_keys:
            label_keys["restraint_type_symbol"] = restraint_type_to_symbol(restraint_type)

        self.state.set_distance_style(start_selector=start_atom,
                                 end_selector=end_atom,
                                 distance_params=distance_params,
                                 sub_style=sub_style,
                                 **label_keys)

        self.state.set_macromolecule_style(
                                selector=self.macromolecule_selector,
                                opacity_params=macromolecule_opacity_params)


    def apply_sub_styles(self, sub_style_func, **sub_style_func_kwargs):
        """
        Apply a sub_style function to annotate all restraints in the restraint_df attribute with their matching sub_style.

        See :ref:`ihm_vis.style.sub_style_modes` for more details on sub_style.

        Parameters
        ----------
        sub_style_func : str or callable
            Either the name of a built-in sub_style mode or a user-provided function
            that accepts (df, **kwargs) and returns a modified DataFrame with a new column defining the sub_style for each row.
        **sub_style_func_kwargs
            Arguments forwarded to the sub_style function.

        Returns
        -------
        df : pandas.DataFrame
            The restraint DataFrame updated with new sub_style assignments.

        Raises
        ------
        ValueError
            If `sub_style_func` (when string) is not found among built_ins.

        """

        # Avoid circular imports
        from ihm_vis.style.sub_style_modes import BUILTIN_SUB_STYLE_FUNCS

        if isinstance(sub_style_func, str):
            sub_style_func = BUILTIN_SUB_STYLE_FUNCS[sub_style_func]

        if isinstance(sub_style_func, str):
            _sub_style_func = BUILTIN_SUB_STYLE_FUNCS.get(sub_style_func, None)

            if _sub_style_func is None:
                raise ValueError(f"The requested built-in sub_style_func ({sub_style_func}) could not be found. See ihm_vis.style.sub_style_modes")

        else:
            _sub_style_func = sub_style_func

        if sub_style_func_kwargs is None:
            sub_style_func_kwargs = {}

        self.restraint_df = _sub_style_func(self.restraint_df, **sub_style_func_kwargs)

        return self.restraint_df


    def set_all_restraint_styles(self, sub_style_col: str="sub_style", **kwargs):
        """
        Loop through all restraints and apply :ref:`ihm_vis.IHM_builder.set_single_restraint_style` to each.

        If no sub_style column is present, all restraints will recieve "default" sub_style.

        Parameters
        ----------
        sub_style_col : str, optional
            Column in `self.restraint_df` that indicates which sub_style to use
            for each restraint (default 'sub_style').
        **kwargs : dict
            Passed through to every call of `set_single_restraint_style`. This can be used to change style for all restraints (as opposed to a sub_set, which is achieve with the sub_styles).
        """

        if not sub_style_col in self.restraint_df.columns:
            restraint_df[sub_style_col] = "default"

        for idx, row in self.restraint_df.iterrows():

            restraint_info = {
                "start_asym_id"  : row["asym_id_1"],
                "start_seq_id"   : row["seq_id_1"], "start_atom_id"  : row["atom_id_1"],

                "end_asym_id"    : row["asym_id_2"],
                "end_seq_id"     : row["seq_id_2"],
                "end_atom_id"    : row["atom_id_2"],

                "distance"       : row["distance_threshold"],
                "restraint_type" : row["restraint_type"],

                "sub_style"      : row[sub_style_col],
            }

            self.set_single_restraint_style(**restraint_info, **kwargs)

    ###############################################################################################
    # Writing mvsj output files
    ###########################

    @classmethod
    def visualize_component(cls, structure: mvs.builder.Structure, component: ComponentStyle) -> mvs.builder.Representation:
        """
        Call underlying MolViewSpec to visualize a styled component.

        Parameters
        ----------
        structure : mvs.builder.Structure
            MolViewSpec structure object.
        component : ihm_vis.style.ComponentStyle
            A selector with representation/color/opacity parameters.

        Returns
        -------
        rep : mvs.builder.Representation
            MolViewSpec representation.
        """

        rep = structure.component(selector=component.selector).representation(**component.representation_params)
        
        if component.color_params:
            rep.color(**component.color_params)

        if component.opacity_params:
            rep.opacity(**component.opacity_params)

        return rep

    @classmethod
    def visualize_distance(cls, structure: mvs.builder.Structure, distance: DistanceStyle) -> mvs.builder.Representation:
        """
        Call underlying MolViewSpec to visualize a distance primitive.

        Parameters
        ----------
        structure : mvs.builder.Structure
            MolViewSpec structure object.
        distance : ihm_vis.style.DistanceStyle
            A distance state containing selectors, params, and label_keys.

        Returns
        -------
        prim : mvs.builder.Representation
            MolViewSpec distance primitive representation
        """

        if "label_template" in distance.distance_params:
            distance.distance_params["label_template"] = distance.distance_params["label_template"].format(**distance.label_keys)

        prim = structure.primitives().distance(start=distance.start_selector, end=distance.end_selector,
                                               **distance.distance_params)

        return prim


    def to_mvsj(self, destination: str|Path, title: Optional[str]="", **kwargs) -> Path:
        """
        Render the current state to a MolViewSpec JSON (.mvsj or .mvsx) file.

        If the source structure was a local file, the mvsj file will be modified so that that local cif file can be bundled with the mvsj into a mvsx archive. This ensures Mol* will have access to the underlying structure.

        Parameters
        ----------
        destination: str or Path
            output file name. Will update extension to 'mvsx' if local file source.
        title : str, optional
            Title metadata to embed in the scene.
        **kwargs
            Additional arguments forwarded to `mvs_builder.save_state`.

        Returns
        -------
        out_path : pathlib.Path
            Path to the generated .mvsj (or .mvsx, if local) file.

        """

        destination = Path(destination)
        if not destination.suffix:
            destination = destination.with_suffix(".mvsj")

        mvs_builder = mvs.create_builder()
        structure = mvs_builder.download(url=self.url).parse(format=self.format).assembly_structure(model_index=self.structure_index)

        # Macromolecule
        self.visualize_component(structure, self.state.macromolecule)

        # Components
        for component in self.state.components.values():
            self.visualize_component(structure, component)

        # Distances
        for distance in self.state.distances.values():
            self.visualize_distance(structure, distance)

        # Write
        if self.source_type == "file":
            archive = destination.with_suffix(".mvsx")

            with self.local_file.serve():
                mvs_builder.save_state(destination=destination, title=title, **kwargs)
                mvsj_to_mvsx(destination, archive, download_external=True)

            return archive

        else:
            mvs_builder.save_state(destination=destination, title=title, **kwargs)

            return destination


