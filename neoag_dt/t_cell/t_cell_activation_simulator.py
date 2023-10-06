""" This module contains a class for determining activated T cells in response
to a vaccine.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List, Optional, Sequence

import numpy as np
import pyllars.collection_utils as collection_utils

from neoag_dt.cell.peptide import Peptide
from neoag_dt.t_cell.t_cell import TCell
from neoag_dt.peptide_score_cache import PeptideScoreCache
from neoag_dt.peptide_vaccine import PeptideVaccine


_DEFAULT_NAME = "TCellActivationSimulator"


class TCellActivationSimulator(object):
    """ A class to determine activated T cells based on a vaccine

    Attributes
    ----------
    dfs_cache : neoag_dt.PeptideScoreCache
        A cache containing the distance-from-self value for all peptides

    name : str
        A name for this object
    """
    def __init__(self,
            dfs_cache:PeptideScoreCache,
            name:str=_DEFAULT_NAME) -> "TCellActivationSimulator":

        self.dfs_cache = dfs_cache
        self.name = name

    def log(self, msg:str, level:int=logging.INFO):    
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def _check_activation(self, peptide:Peptide) -> Optional[TCell]:
        """ Check if `peptide` activates any T cells

        Parameters
        ----------
        peptide : neoag_dt.Peptide
            The peptide

        Returns
        -------
        activated_tcell : typing.Optional[neoag_dt.TCell]
            The activated T cell, if `peptide` led to activation. If it did not,
            then `None` is returned.
        """
    
        dfs_score = self.dfs_cache.get_score(peptide)
        t_cell_is_activated = np.random.binomial(n=1, p=dfs_score)

        activated_t_cell = None
        if t_cell_is_activated:
            known_targets = [peptide]
            activated_t_cell = TCell(known_targets=known_targets)

        return activated_t_cell

    def simulate(self, vaccine:PeptideVaccine) -> List[TCell]:
        """ Based on the contents of `vaccine`, determine the set of activated
        T cells

        Parameters
        ----------
        vaccine : neoag_dt.PeptideVaccine
            The vaccine

        Returns
        -------
        activated_t_cells : typing.Sequence[neoag_dt.TCell]
            The set of T cells activated by `vaccine`
        """
        activated_t_cells = collection_utils.remove_nones([
            self._check_activation(peptide) for peptide in vaccine.peptides
        ])

        return activated_t_cells

