"""
"""

import molviewspec as mvs
from typing import Optional, Tuple, Dict, Callable, Any

from ihm_vis.style import DEFAULT, apply_style_defaults
from ihm_vis.sub_style_modes import BUILTIN_SUBSTYLE_FUNCS


class IHM_Builder:

    restraint_df_columns  = ["entity_id_1", "asym_id_1", "seq_id_1", "comp_id_1", "atom_id_1",
                             "entity_id_2", "asym_id_2", "seq_id_2", "comp_id_2", "atom_id_2",
                             "model_granularity", "distance_threshold", "restraint_type",
                            ]

    def __init__(self, source: str, mvs_builder: Optional[mvs.Builder]=None, structure_index: int=0, format: str="mmcif"):
        """
        """

        # Set up MolViewSpec underlying
        # builder and basic environment
        if mvs_builder is not None:
            self.mvs_builder = mvs_builder
        else:
            self.mvs_builder = mvs.create_builder()

        self.structure = self.mvs_builder.download(url=url).parse(format=format).assembly_structure()

        # Parse cif file for restraint information
        self.cif = read_cif(source, structure_index)

        # Setup for restraint_df
        self._restraint_df_init = False


    ###############################################################################################
    # General Utility functions
    ###########################

    @classmethod
    def read_cif(cls, source: str, structure_index: int=0):  # TODO: mmcif type
        """
        """
        ...

    @classmethod
    def get_atom_coordinates(cls, atom_id, comp_id, entity_id, asym_id, seq_id) -> Optional[Tuple[float, float, float]]:
        """
        Retrieve the Cartesian coordinates (x, y, z) of an atom based on its identifiers.

        Parameters
        ----------
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
        Extract cross-link restraint data from a list of containers and return it as a pandas DataFrame.

        Returns
        -------
        df : pd.DataFrame
            A pandas DataFrame where each row represents a cross-link restraint, ith columns including:
            'entity_id_1', 'asym_id_1', 'seq_id_1', 'comp_id_1', 'atom_id_1',
            'entity_id_2', 'asym_id_2', 'seq_id_2', 'comp_id_2', 'atom_id_2',
            'model_granularity', 'distance_threshold', and 'restraint_type'.
        """

        # Extract cross-link restraint data from the container
        cross_link_data = self.cif.getObj('ihm_cross_link_restraint')
        if cross_link_data is None:
            return pd.DataFrame(columns=self.restraint_df_columns)

        # Iterate through each row of cross-link data
        restraints = []
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
            "restraint_type": cross_link_data.getValue("restraint_type", i),
            }
            restraints.append(data)
    
        # Convert the list of restraints to a pandas DataFrame.
        restraint_df = pd.DataFrame(restraints)
        # If atom ids are not specified default to carbon alpha (CA)
        restraint_df['atom_id_1'] = restraint_df['atom_id_1'].str.replace('.', 'CA', regex=False)
        restraint_df['atom_id_2'] = restraint_df['atom_id_2'].str.replace('.', 'CA', regex=False)
        # Lists to store retrieved atom coordinates
        atom_id_1_coords = []
        atom_id_2_coords = []
        # Iterate through cross link in dataframe
        for index, row in self.restraint_df.iterrows():
            # Get coordinates for each atom in crosslink
            atom_id_1_coords.append(get_atom_coordinates(container, row['atom_id_1'], row['comp_id_1'], row['entity_id_1'], row['asym_id_1'], row['seq_id_1']))
            atom_id_2_coords.append(get_atom_coordinates(container, row['atom_id_2'], row['comp_id_2'], row['entity_id_2'], row['asym_id_2'], row['seq_id_2']))
        # Add coordinates to dataframe
        restraint_df['atom_id_1_coords'] = atom_id_1_coords
        restraint_df['atom_id_2_coords'] = atom_id_2_coords

        self.restraint_df = pd.concat((self.restraint_df, restraint_df)).drop_duplicates()
        self._restraint_df_init = True
        return restraint_df

    
    # Master method to parse all restraint types #
    ##############################################
    @property
    def restraint_df(self) -> pd.DataFrame:
        if self._restraint_df_init:
            return self.restraint_df

        else:
            # Call each restraint type
            self.get_cross_links()


        self._restraint_df_init = True
        return self.restraint_df


    ##############################################################################################
    # Filtering restraints
    ######################

    def filter_restraints(filter_func: str|Callable[[IHM_Builder], pd.DataFrame], **kwargs) -> pd.DataFrame:
        """
        """
        if isinstance(filter_func, str):
            _filter_func = BUILTIN_FILTER_FUNCS.get(filter_func, None)

            if _filter_func is None:
                raise ValueError(f"The requested built-in filter_func ({filter_func}) could not be found. See ihm_vis.filters")

        else:
            _filter_func = filter_func

        self.restraint_df = _filter_func(self, **kwargs)
        return self.restraint_df



    ###############################################################################################
    # Visualization Functions
    #########################

    @apply_style_defaults
    def visualize_macromolecule(self, 
                                representation_params: Optional[Dict[str, str]=DEFAULT, 
                                color_params: Optional[Dict[str, str]]=DEFAULT, 
                                opacity_params: Optional[Dict[str, str]]=DEFAULT):
		"""
		Visualize the macromolecule structure

		Parameters
		----------
			representation_params: representation parameters passed to MolViewSpec
			color_params: color parameters passed to MolViewSpec
            opacity_params: opacity parameters passed to MolViewSpec

		Raises
		------
            TypeError: Parameter passed to MolViewSpec is not supported
		"""
        self.structure.component(selector="polymer").representation(**representation_params).color(**color_params).opacity(**opacity_params)


	@apply_style_defaults
	def visualize_restraint(self,

							start_asym_id: str|int, 
							start_seq_id: int, 

							end_asym_id: str|int,
							end_seq_id: int, 

							distance: float,
							restraint_type: str, 

							start_atom_id: str="CA", end_atom_id: str="CA",

							representation_params: Optional[Dict[str, str]]=DEFAULT,
							color_params: Optional[Dict[str, str]]=DEFAULT,
							distance_params: Optional[Dict[str, str]]=DEFAULT,
							tube_params: Optional[Dict[str, str]]=DEFAULT,

							sub_style="default",

							focus: Optional[bool]=False):

		"""
		Visualize a restraint

		Parameters
		----------

			start_asym_id: str | int, asym_id of the starting residue
			start_seq_id: int, seq_id of the starting residue

			end_asym_id: str | int, asym_id of the ending residue
			end_seq_id: int, seq_id of the ending residue

			distance: float, the restraint distance
			restraint_type: str, the type/operator of the restraint, one of {list(RESTRAINT_TYPE_TO_SYMBOL.keys())}

			representation: str="ball_and_stick", representation of the two residues of the restraint
			residue_color: str="red", color of the restraint residues
			start_atom_id: str="CA", atom_id of the starting residue
			end_atom_id: str="CA", atom_id of the ending residue

			radius: float=0.1, radius of the line drawn between restraint residues
			line_color=None, color for the line drawn between restraint residues, defaults to the residues' color
			dash_length: float=0.1, dash length of the line drawn between restraint residues

			label_template="Solved Distance: {{{{distance}}}}, Restraint Distance: the text displayed on the line drawn between restraint residues
			label_color=None, color of the text displayed on the line drawn between residues, defaults to the residues' color

			focus: bool=True, whether to focus the carmera to this restraint

		Raises
		------
            TypeError: Parameter passed to MolViewSpec is not supported
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

		start_component = self.structure.component(selector=start_residue)
		start_component.representation(**representation_params).color(**color_params).opacity(**opacity_params)

		end_component = self.structure.component(selector=end_residue)
		end_component.representation(**representation_params).color(**color_params).opacity(**opacity_params)

		if distance_params is not None:
			if "tooltip" in distance_params:
				distance_params["tooltip"] = distance_params["tooltip"].format(restraint_type_symbol=restraint_type_to_symbol(restraint_type), distance=distance)

			res = self.structure.primitives().distance(
					start=start_atom,
					end=end_atom,
					**distance_params)

		if tube_params is not None:
			if "tooltip" in tube_params:
				tube_params["tooltip"] = tube_params["tooltip"].format(restraint_type_symbol=restraint_type_to_symbol(restraint_type), distance=distance)
			res = self.structure.primitives().tube(
					start=start_atom,
					end=end_atom,
					**tube_params)

		if focus:
			res.focus()


    def visualize_restraints(self, sub_style_func: Optional[str|Callable[[IHM_Builder], pd.Series]]="default", sub_style_func_kwargs: Optional[Dict[Any, Any]]=None, **kwargs):
		"""
		Visualize all restraints.

        This function loops over all restraints in the ihm_builder.restraint_df. If no specific restraints tpyes have been parsed yet,
        such as by explicityly calling ihm_builder.get_cross_links(), then all restraint types will be parsed before visualization

		Parameters
		----------
            sub_style_func: Optional[str|Callable[[IHM_Builder], pd.DataFrame]], a user-defined function that given the builder
                             will return a pandas series of sub_styles for each row of the restraint_df
                             built-in functions for common operations provided in ihm_vis.sub_style_modes.

            sub_style_func_kwargs: Optional[Dict[Any, Any]], additional arguments to pass to sub_style_func

            **kwargs: Any arguments to pass to IHM_Builder.visualize_restraint for every row

		Raises
		------
            TypeError: Parameter passed to MolViewSpec is not supported
		"""

        if isinstance(sub_style_func, str):
            sub_style_func = BUILTIN_SUBSTYLE_FUNCS[sub_style_func]

        if isinstance(sub_style_func, str):
            _sub_style_func = BUILTIN_SUB_STYLE_FUNCS.get(sub_style_func, None)

            if _sub_style_func is None:
                raise ValueError(f"The requested built-in sub_style_func ({sub_style_func}) could not be found. See ihm_vis.sub_style_modes")

        else:
            _sub_style_func = sub_style_func

        if sub_style_func_kwargs is None:
            sub_style_func_kwargs = {}

        sub_styles = _sub_style_func(self, **sub_style_func_kwargs)

        for idx, row in self.restraint_df.iterrows():

            restraint_info = {
                "start_asym_id"  : row["asym_id_1"],
                "start_seq_id"   : row["seq_id_1"], "start_atom_id"  : row["atom_id_1"],

                "end_asym_id"    : row["asym_id_2"],
                "end_seq_id"     : row["seq_id_2"],
                "end_atom_id"    : row["atom_id_2"],

                "distance"       : row["distance_threshold"],
                "restraint_type" : row["restraint_type"],

                "sub_style"      : sub_styles.loc[idx],
            }

            self.visualize_restraint(structure, **restraint_info, **kwargs)


