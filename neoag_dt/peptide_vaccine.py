""" This module contains a class for representing vaccines.

Specifically, it considers a vaccine as a set of peptides. It does not consider
spacers or other factors relevant for actual vaccines.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Sequence
from neoag_dt.cell.peptide import Peptide


_DEFAULT_NAME = "PeptideVaccine"


class PeptideVaccine(object):
    """ A vaccine which considers only peptide sequences
    """
    def __init__(self,
            peptides:Sequence[Peptide],
            name:str=None) -> "PeptideVaccine":

        self.name = name
        self._peptides = peptides

    def log(self, msg:str, level:int=logging.INFO):    
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    @property
    def peptides(self):
        return self._peptides

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.name

