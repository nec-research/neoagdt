""" This module contains code necessary for peptide vaccine design based on
our bipartite (peptides - cells) formulation.
"""
import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import tqdm

import mip

from neoag_dt.optimization import optimization_utils

from typing import Mapping, Sequence

_DEFAULT_NAME = "BipartiteVaccineDesignModel"
_DEFAULT_MAX_SOLVING_TIME = 3600  # an hour


class BipartiteVaccineDesignModel(object):
    """ A class which implements the max flow problem for
    optimal vaccine design.

    Attributes
    ----------

    config : Mapping
        The optimization configuration.

    peptides : Sequence[str]
        The candidates.

    df_cells : pd.DataFrame
        The cancer digital twin.

    df_scores : pd.DataFrame
        The simulated response likelihood.

    peptide_weights_map : Mapping[str, float]
        The weight/cost of each vaccine element.

    vaccine_elements : str
        The vaccine elements.

    seed : int
        Random seed.

    name : str
        Object name.
    """
    def __init__(
            self,
            config: Mapping,
            peptides: Sequence[str],
            df_cells: pd.DataFrame,
            df_scores: pd.DataFrame,
            peptide_weights_map: Mapping[str, float],
            vaccine_elements: str,
            seed: int,
            name: str = _DEFAULT_NAME):

        self.config = config
        self.peptides = peptides
        self.df_cells = df_cells
        self.df_scores = df_scores
        self.peptide_weights_map = peptide_weights_map
        self.vaccine_elements = vaccine_elements
        self.seed = seed
        self.name = name

        self._initialize()

        self._status = None

    def _initialize(self):

        # peptide maps and counts
        peptide_maps = optimization_utils.get_index_and_reverse_map(
            self.peptides
        )
        self.peptide_index_map, self.peptide_reverse_index_map = peptide_maps
        self.num_peptides = len(self.peptide_index_map)

        # cell maps and counts
        cell_maps = optimization_utils.get_index_and_reverse_map(
            self.df_cells['cell_ids'].unique()  # only consider APCs
        )
        self.cell_index_map, self.cell_reverse_index_map = cell_maps
        self.num_cells = len(self.cell_index_map)

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def _add_variables(self):
        # add the decision variables for each peptide
        self.x_peptide = [
            self.m.add_var(var_type=mip.BINARY, name="x_peptide_{}".format(i))
            for i in range(self.num_peptides)
        ]
        # add the continuous variables for each cell
        self.x_cell = [
            self.m.add_var(var_type=mip.CONTINUOUS, lb=-np.inf, name="x_cell_{}".format(k))
            for k in range(self.num_cells)
        ]

        # add the continuous variable we will minimize
        self.z = self.m.add_var(var_type=mip.CONTINUOUS, lb=-np.inf, name="z")

    ###
    # Peptide-allele scores and constraints
    ###
    def _create_peptide_cell_scores_matrix(self):
        self.peptide_cell_scores = optimization_utils.get_peptide_cell_scores_matrix(
            self.df_scores, self.peptide_index_map, self.cell_index_map
        )

    def _add_cell_score_constraint(self, cell):
        cell_index = self.cell_index_map[cell]
        # the logP(R=- | V, c_j) of no response for a cell c_j and
        # a given vaccine V is the sum of the scores logP(R=- | v_i, c_j)
        # of no response for each vaccine element v_i
        cell_log_p_no_response = mip.xsum(
            self.peptide_cell_scores[i, cell_index] * self.x_peptide[i]
            for i in range(self.num_peptides)
        )

        constraint_name = "cell_log_p_no_resopnse_{}".format(cell)
        constraint = (cell_log_p_no_response == self.x_cell[cell_index])

        self.m += (constraint, constraint_name)

    def _add_cell_score_constraints(self):
        cells = self.cell_index_map.keys()

        for cell in tqdm.tqdm(cells):
            self._add_cell_score_constraint(cell)

    ###
    # Other constraints
    ###
    def _add_fixed_peptide_constraints(self):
        """ Constraint which ensures that a set of peptides is included
         in the vaccine design.
         """
        fixed_peptides = self.config.get('fixed_peptides')
        if fixed_peptides is None:
            return

        for peptide in fixed_peptides:
            peptide_index = self.peptide_index_map[peptide]
            var_name = "x_peptide_{}".format(peptide_index)
            variable = self.m.var_by_name(var_name)

            if variable is None:
                msg = "Could not find fixed variable: '{}'".format(peptide)
                raise ValueError(msg)

            constraint_name = "fixed_x_peptide_{}".format(peptide)
            constraint = (variable == 1)
            self.m += (constraint, constraint_name)

    def _add_budget_constraint(self):
        self.peptide_weights = [
            self.peptide_weights_map[self.peptide_reverse_index_map[i]]
            for i in range(self.num_peptides)
        ]

        total_cost = mip.xsum(
            self.peptide_weights[i] * self.x_peptide[i]
            for i in range(self.num_peptides)
        )

        self.m += (total_cost <= self.config['budget'], 'budget')

    def _add_minimax_constraint(self):
        # the constraint that z must be greater than all x_cell[j] values
        for j in range(self.num_cells):
            self.m += (self.x_cell[j] <= self.z, 'cell_to_z_{}'.format(j))

    def _add_minsum_constraint(self):
        log_p_no_response = mip.xsum(
            self.x_cell[j] for j in range(self.num_cells)
        )

        constraint_name = "log_p_no_resopnse_{}"
        constraint = (log_p_no_response == self.z)
        self.m += (constraint, constraint_name)

    def _add_constraints(self):
        """ Add the ILP constraints.
        """
        self._add_cell_score_constraints()
        self._add_budget_constraint()
        self._add_fixed_peptide_constraints()

        if self.config['criterion'] == 'MinMax':
            self._add_minimax_constraint()
        elif self.config['criterion'] == 'MinSum':
            self._add_minsum_constraint()
        else:
            raise NotImplementedError(f"{self.config['criterion']} is not implemented")

    def build_model(self):
        """ Build the bipartite graph using the
        `mip` library.
        """
        msg = "creating peptide cell scores matrix"
        self.log(msg)
        self._create_peptide_cell_scores_matrix()

        msg = "creating model"
        self.log(msg)
        self.m = mip.Model(self.name)
        self.m.seed = self.seed

        msg = "adding variables to model"
        self.log(msg)
        self._add_variables()

        msg = "adding constraints to model"
        self.log(msg)
        self._add_constraints()

        # and set the optimization problem
        self.m.objective = mip.minimize(self.z)

    def optimize(self, max_solving_time=_DEFAULT_MAX_SOLVING_TIME) -> "BipartiteVaccineDesignModel":
        """ Run optimization.
        """
        msg = "optimizing the peptide selection"
        self.log(msg)
        self._status = self.m.optimize(max_seconds=max_solving_time)
        return self

    def get_selected_peptides(self) -> pd.DataFrame:
        """ Get the best vaccine elements.
        """
        selected_peptides = []

        for i in range(self.num_peptides):
            var_name = "x_peptide_{}".format(i)
            x_peptide_i = self.m.var_by_name(var_name)
            selected = (x_peptide_i.x > 0.99)

            if selected:
                peptide = self.peptide_reverse_index_map[i]
                selected_peptides.append(peptide)

        repetition = self.df_cells.iloc[0]['repetition']
        simulation_name = self.df_cells.iloc[0]['simulation_name']

        df_selected_peptides = pd.DataFrame()
        df_selected_peptides['peptide'] = selected_peptides
        df_selected_peptides['repetition'] = repetition
        df_selected_peptides['simulation_name'] = simulation_name

        return df_selected_peptides

    def check_optimization(self) -> None:
        """ Check the optimization has been correctly
        solved.
        """
        self.m.check_optimization_results()

        if self._status == mip.OptimizationStatus.OPTIMAL:
            logger.info('optimal solution cost {} found'.format(
                self.m.objective_value
            ))
        elif self._status == mip.OptimizationStatus.FEASIBLE:
            logger.info('sol.cost {} found, best possible: {}'.format(
                self.m.objective_value, self.m.objective_bound
            ))
        elif self._status == mip.OptimizationStatus.NO_SOLUTION_FOUND:
            logger.info('no feasible solution found, lower bound is: {}'.format(
                self.m.objective_bound
            ))
