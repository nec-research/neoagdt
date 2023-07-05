""" This module contains a class for representing a cancer cell.
"""
import logging
logger = logging.getLogger(__name__)

from neoag_dt.cell.genetic_simulator import GeneticSimulationResult
from neoag_dt.cell.mhc import MHC
from neoag_dt.cell.peptide import Peptide
from neoag_dt.cell.somatic_variant import SomaticVariant
from neoag_dt.cell.peptide_mhc_binding_simulator import pMHC

from typing import Mapping, Optional, Sequence, Set


class Cell(object):
    """ A cancer cell.

    **N.B.** All of the attributes are optional. This is largely to reflect
    that different simulations may not require all attributes.

    Attributes
    ----------
    genetic_simulation_results : Mapping[SomaticVariant, GeneticSimulationResult]
        The results of the simulations for each candidate somatic variant. This
        includes whether the variant was present in the DNA and how many protein
        molecules containing the variant are estimated to be present.

    selected_peptides : Sequence[Peptide]
        The peptides which were selected

    hla_counts : Mapping[MHC, int]
        The number of each HLA molecule present inside the cell

    bound_peptides : Sequence[pMHC]
        The peptide:MHC complexes which are bound together in the inside of the
        cell. Not all of the bound peptides are necessarily presented on the
        cell surface.

    presented_peptides : Sequence[pMHC]
        The peptide:MHC complexes presented on the cell surface

    name : str
        A name for this cell
    """

    def __init__(
            self,
            genetic_simulation_results: Optional[Mapping[SomaticVariant, GeneticSimulationResult]] = None,
            selected_peptides: Optional[Sequence[Peptide]] = None,
            hla_counts: Optional[Mapping[MHC, int]] = None,
            bound_peptide_hlas: Optional[Sequence[pMHC]] = None,
            presented_peptide_hlas: Optional[Sequence[pMHC]] = None,
            name: Optional[str] = None) -> "Cell":

        self.name = name
        self.genetic_simulation_results = genetic_simulation_results
        self.selected_peptides = selected_peptides
        self.hla_counts = hla_counts
        self.bound_peptide_hlas = bound_peptide_hlas
        self.presented_peptide_hlas = presented_peptide_hlas

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    @property
    def presented_peptides(self) -> Set[Peptide]:
        presented_peptides = {
            p.peptide for p in self.presented_peptide_hlas
        }
        return presented_peptides

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.name
