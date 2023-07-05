""" This module contains a class for vaccine evaluation.

N.B. Likelihood of no response for a cell c_j, given a vaccine V:
    P(R =  - | V, c_j) = \prod_{i=1}^{N} P(R =  - | v_i, c_j)
    log P(R =  - | V, c_j) = \sum_{i=1}^{N} log P(R =  - | v_i, c_j)
"""
import logging
logger = logging.getLogger(__name__)

import numpy as np
import pyllars.validation_utils as validation_utils
import tqdm
import pandas as pd

from neoag_dt.evaluation.cell import Cell
from neoag_dt.evaluation.population import Population
from neoag_dt.evaluation.vaccine import Vaccine

_DEFAULT_NAME = "VaccineEvaluator"


class VaccineEvaluator:
    """ A class for a vaccine evaluator.

    Attributes
    ----------
    population : Population
        The cell population.

    vaccine : Vaccine
        A vaccine object (set of vaccine elements).

    df_scores : pd.DataFrame
        Response likelihood estimated by the optimization.

    name : str
        Object name.
    """
    def __init__(
            self,
            population: Population,
            vaccine: Vaccine,
            df_scores: pd.DataFrame,
            name: str = _DEFAULT_NAME) -> "VaccineEvaluator":

        self.population = population
        self.vaccine = vaccine
        self.df_scores = df_scores
        self.name = name

        self._create_scores_maps()

        self._p_no_response = dict()
        self._log_p_no_response = dict()
        self._p_response = dict()
        self._log_p_response = dict()

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger
        """
        msg = "[{}: {}] {}".format(self.name, self.vaccine.name, msg)
        logger.log(level, msg)
        
    def get_p_response_values(self):
        """ Get response probabilities.
        """
        validation_utils.check_is_fitted(self, '_p_response')
        ret = np.array(list(self._p_response.values()))
        return ret

    def get_log_p_response_values(self):
        """ Get log-probabilities of response.
        """
        validation_utils.check_is_fitted(self, '_log_p_response')
        ret = np.array(list(self._log_p_response.values()))
        return ret

    def _evaluate_cell(self, cell: Cell) -> None:
        """ Computes log P(R =  - | V, c_j) and stores response scores in
        related dictionaries."""
        log_p_no_response_cell = self._get_cell_response(cell.name)
        self._store_cell_scores(cell, log_p_no_response_cell)

    def _get_cell_response(self, cell_id: int) -> float:
        """ Computes the probability of no response for a given cell.
        N.B.
            log P(R =  - | V, c_j) = \sum_{i=1}^{N} log P(R =  - | v_i, c_j)
        """
        scores = [
            self._scores_maps[cell_id][v]
            for v in self.vaccine
        ]

        # we assume vaccine elements are independent, so the log likelihood that all of
        # them fail is just the sum
        log_p_no_response = np.sum(scores)

        return log_p_no_response

    def _store_cell_scores(
            self,
            cell: Cell,
            log_p_no_response_cell: float) -> None:
        """ Stores response scores for a cell in related dictionaries."""
        p_no_response_cell = np.exp(log_p_no_response_cell)
        p_response_cell = 1 - p_no_response_cell
        log_p_response_cell = np.log(p_response_cell)

        self._p_no_response[cell.name] = p_no_response_cell
        self._log_p_no_response[cell.name] = log_p_no_response_cell
        self._p_response[cell.name] = p_response_cell
        self._log_p_response[cell.name] = log_p_response_cell

    def _create_scores_maps(self, progress_bar: bool = True) -> None:
        """ Create a list of score maps.
        A score map is created for each cell in the simulation.
        Each score map is a dictionary which maps the vaccine elements
        to the `log_p_no_response` coming from the scores file.
        """
        self._scores_maps = []

        msg = "creating scores maps"
        self.log(msg)
        it = self.population.cells
        if progress_bar:
            it = tqdm.tqdm(it)

        for c in it:
            df_cell = self.df_scores[self.df_scores['cell_id'] == c.name]
            score_map = dict(zip(df_cell['vaccine_element'], df_cell['log_p_no_response']))
            self._scores_maps.append(score_map)

    def evaluate(self, progress_bar: bool = True) -> "VaccineEvaluator":
        """ Evaluates the response likelihood for each cell in a population."""
        msg = "evaluating cells"
        self.log(msg)
        it = self.population.cells
        if progress_bar:
            it = tqdm.tqdm(it)

        for cell in it:
            self._evaluate_cell(cell)

        return self
