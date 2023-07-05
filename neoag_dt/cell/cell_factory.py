""" This module contains a class for simulating cancer cells.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List, Mapping, Optional, Sequence

import numpy as np
import pyllars.collection_utils as collection_utils
import tqdm

from neoag_dt.cell.genetic_simulator import GeneticSimulator
from neoag_dt.cell.mhc import MHC
from neoag_dt.cell.peptide import Peptide
from neoag_dt.cell.somatic_variant import SomaticVariant
from neoag_dt.cell.peptide_mhc_binding_simulator import PeptideMHCBindingSimulator
from neoag_dt.cell.cell import Cell
from neoag_dt.cell.allele_score_cache import AlleleScoreCache
from neoag_dt.cell.protein_cleaver import ProteinCleaver
from neoag_dt.cell.peptide_mhc_binding_simulator import pMHC
import neoag_dt.dt_utils as dt_utils

_DEFAULT_NAME = "CellFactory"


class CellFactory(object):
    """ A class to simulate all of the internal processes within a cancer cell
    to generate the final set of peptide:MHC complexes presented on the outside
    of the cell surface.
    
    Attributes
    ----------
    genetic_simulator : neoag_dt.GeneticSimulator

    protein_cleaver : neoag_dt.ProteinCleaver

    binding_simulator : neoag_dt.PeptideMHCBindingSimulator

    presentation_scores : neoag_dt.AlleleScoreCache
    
    expression_pseudocount : float
        A pseudocount to add to all HLA expression values.These can be thought
        of as a kind of hyperprior on the parameters of the gamma distribution
        for sampling.

        Practically, these are used to avoid passing `0`s to the sampling
        functions.

        It is not straightforward to express the Jeffreys prior for the gamma
        distribution, especially if the observed expression is close to `0`. So
        we just use `1` as the default.

        Please see this discussion for more details on the Jeffreys prior for
        a gamma distribution: https://stats.stackexchange.com/questions/197173

    name : str
        A name for this object
    """
    def __init__(
            self,
            genetic_simulator: GeneticSimulator,
            protein_cleaver: ProteinCleaver,
            binding_simulator: PeptideMHCBindingSimulator,
            presentation_scores: AlleleScoreCache,
            expression_pseudocount: int = 1,
            name=_DEFAULT_NAME) -> "CellFactory":

        self.name = name
        self.genetic_simulator = genetic_simulator
        self.protein_cleaver = protein_cleaver
        self.binding_simulator = binding_simulator
        self.presentation_scores = presentation_scores
        self.expression_pseudocount = expression_pseudocount

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def _get_selected_peptides(self, variant, simulation_result, hla_alleles):
        num_proteins = simulation_result.protein_count
            
        selected_peptides = self.protein_cleaver.cleave(
            variant.peptides, hla_alleles, num_proteins
        )

        return selected_peptides

    def _get_all_selected_peptides(
            self,
            hla_alleles: Sequence[MHC],
            progress_bar: bool) -> Sequence[Peptide]:

        it = self.genetic_simulation_results.items()
        if progress_bar:
            it = tqdm.tqdm(it)

        all_selected_peptides = collection_utils.flatten_lists([
            self._get_selected_peptides(v, r, hla_alleles) for v, r in it
        ])

        return all_selected_peptides

    def _get_hla_counts(self, hla_alleles: Sequence[MHC]) -> Mapping[MHC, int]:

        hla_counts = {
            hla: dt_utils.sample_gamma_poisson(
                hla.protein.expression_mean + self.expression_pseudocount,
                hla.protein.expression_var
            ) for hla in hla_alleles
        }

        return hla_counts

    def _get_presented_peptide_hla_complexes(
            self,
            bound_peptide_hla_complexes: Sequence[pMHC],
            progress_bar: bool = True) -> List[pMHC]:

        # for each bound complex
        #   check if the complex is presented
        presentation_likelihoods = np.array([
            self.presentation_scores.get_score(hla, peptide)
                for peptide, hla in bound_peptide_hla_complexes
        ])

        presented = np.random.binomial(n=1, p=presentation_likelihoods)
        presented = np.where(presented)[0]

        bound_peptide_hla_complexes = np.array(bound_peptide_hla_complexes)
        presented_peptides = bound_peptide_hla_complexes[presented]
        presented_peptides = presented_peptides.tolist()

        #TODO: this is inefficent
        # convert the presented peptides back to pMHC objects
        presented_peptides = [
            pMHC(
                peptide=p[0],
                hla=p[1]
            ) for p in presented_peptides
        ]

        return presented_peptides

    def create_cell(
            self,
            variant_map: Mapping[str, SomaticVariant],
            hla_alleles: Sequence[MHC],
            name: Optional[str] = None,
            progress_bar: bool = True) -> Cell:

        # simulate the genetic impact of each variant
        self.genetic_simulation_results = self.genetic_simulator.simulate_all(variant_map, progress_bar)
            
        # determine the peptides from each variant
        self.selected_peptides = self._get_all_selected_peptides(hla_alleles, progress_bar)

        # number of HLA molecules
        self.hla_counts = self._get_hla_counts(hla_alleles)

        # binding
        self.bound_peptides = self.binding_simulator.get_bound_peptide_hla_complexes(
            self.selected_peptides, self.hla_counts, progress_bar
        )

        # presentation
        self.presented_peptides = self._get_presented_peptide_hla_complexes(
            self.bound_peptides, progress_bar
        )

        cell = Cell(
            genetic_simulation_results=self.genetic_simulation_results,
            selected_peptides=self.selected_peptides,
            hla_counts=self.hla_counts,
            bound_peptide_hlas=self.bound_peptides,
            presented_peptide_hlas=self.presented_peptides,
            name=name
        )

        return cell
