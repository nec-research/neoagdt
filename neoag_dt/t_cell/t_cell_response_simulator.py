""" This module contains a class for determine whether a set of activated
T cells respond to any of a set of presented peptide:HLA complexes.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List, NamedTuple, Optional, Sequence

import numpy as np
import pyllars.collection_utils as collection_utils

from neoag_dt.cell.cell import Cell
from neoag_dt.t_cell.t_cell import TCell
from neoag_dt.peptide_score_cache import PeptideScoreCache
from neoag_dt.cell.peptide_mhc_binding_simulator import pMHC


_DEFAULT_NAME = "TCellResponseSimulator"


class TCellSimulationResults(NamedTuple):
    peptide_hla: pMHC
    t_cell : TCell
    response_likelihood: float
    response : bool


class TCellResponseSimulator(object):
    """ A class to determine whether T cells respond to pHLAs

    Attributes
    ----------
    t_cell_repertoire : typing.Sequence[TCell]
        The repertoire of activated T cells

    response_cache : neoag_dt.PeptideScoreCache
        A cache containing the likelihood of response for all peptides

    name : str
        A name for this object

    simulation_results : typing.Dict[Cell, typing.Sequence[TCellSimulationResults]]
        A mapping from each cell to all simulation results for that cell. This
        is created as a side effect during the `simulate` method.
    """
    def __init__(self,
            t_cell_repertoire:Sequence[TCell],
            response_cache:PeptideScoreCache,
            name:str=_DEFAULT_NAME) -> "TCellResponseSimulator":

        self._t_cell_repertoire = t_cell_repertoire
        self._response_cache = response_cache
        self.name = name

    def log(self, msg:str, level:int=logging.INFO):    
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def _check_recognition(self, peptide_hla:pMHC) -> List[TCell]:
        """ Check if `pHLA` is recognized by any cells in the T-cell repertoire

        This *does not* check if those T cells actually respond. It only checks
        which among the repertoire recognize the pHLA complex.

        Parameters
        ----------
        pHLA : neoag_dt.pMHC
            The peptide HLA complex

        Returns
        -------
        recognized_tcells : typing.Sequence[neoag_dt.TCell]
            The T cells which recognize the pHLA complex.
        """
        recognized_tcells = [
            t_cell for t_cell in self._t_cell_repertoire
                if t_cell.recognizes_target(peptide_hla)
        ]

        return recognized_tcells

    def _check_response(self, peptide_hla:pMHC, t_cell:TCell) -> TCellSimulationResults:
        """ Check whether `t_cell` responds to `peptide_hla`
        
        Currently, this does not take into account the HLA or the T cell. It
        only checks the likelihood of response based on the peptide sequence,
        and it randomly checks for response.
        """
        peptide = peptide_hla[0]
        response_likelihood = self._response_cache.get_score(peptide)
        response = np.random.binomial(n=1, p=response_likelihood)

        res = TCellSimulationResults(peptide_hla, t_cell, response_likelihood, response)
        return res

    def _simulate_cell(self, cell:Cell) -> List[TCellSimulationResults]:
        """ Check whether any T cell recognizes and responds to `cell` """

        # for each peptide_hla presented on the cell surface, check if any
        # T cell in the T-cell repertoire recognizes that pMHC complex
        all_responses = []


        # for each peptide:HLA complex, check how many T cells recognize it
        for peptide_hla in cell.presented_peptide_hlas:
            recognized_t_cells = self._check_recognition(peptide_hla)

            responses = [
                self._check_response(peptide_hla, t_cell)
                    for t_cell in recognized_t_cells
            ]

            all_responses.append(responses)

        all_responses = collection_utils.flatten_lists(all_responses)
        return all_responses
            

    def simulate(self, cells:Sequence[Cell]) -> List[Cell]:
        """ Based on the T-cell repertoire, determine which cells experience
        an T-cell response

        Presumably, the response indicates that the cell is killed.

        Parameters
        ----------
        cells : typing.Sequence[neoag_dt.Cell]
            All of the simulated cells. Importantly, the `presented_peptides`
            must be set.

        Returns
        -------
        cells_with_response : typing.List[neoag_dt.Cell]
            The cells which elicited a response from at least one T cell

        Notes
        -----
        This method has the side effect of creating a `simulation_results`
        attribute for this object. It is a map from each cell to all simulation
        results for that cell, which include a peptide:MHC complex, a T cell
        which recognized that complex, and the simulation result.
        """
        self.simulation_results = {}
        cells_with_response = []

        for cell in cells:
            cell_simulation_results = self._simulate_cell(cell)
            self.simulation_results[cell] = cell_simulation_results

            cell_response = any(r.response for r in cell_simulation_results)

            if cell_response:
                cells_with_response.append(cell)

        return cells_with_response









    