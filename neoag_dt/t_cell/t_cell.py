""" This module contains a class for representing a T cell.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Sequence

from neoag_dt.cell.peptide import Peptide
from neoag_dt.cell.peptide_mhc_binding_simulator import pMHC


class TCell(object):
    """ A T cell
    """
    def __init__(self,
            known_targets:Sequence[Peptide],
            name:str=None) -> "TCell":

        self._known_targets = known_targets
        self.name = name

        self._known_targets_set = {p for p in known_targets}

    def log(self, msg:str, level:int=logging.INFO):    
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def recognizes_target(self, peptide_hla:pMHC) -> bool:
        """ Check whether `peptide_hla` is recognized for this T cell
        
        Currently, this just checks if the peptide is among the known targets
        for this T cell; the HLA is not considered.
        """
        recognized = peptide_hla.peptide in self._known_targets_set
        return recognized

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.name