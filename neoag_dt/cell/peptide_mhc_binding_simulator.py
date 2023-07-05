""" This module contains a class for simulating peptide:HLA binding.
"""
import logging

logger = logging.getLogger(__name__)

from typing import List, Mapping, Sequence, NamedTuple

import numpy as np
import tqdm
import pyllars.collection_utils as collection_utils

from neoag_dt.cell.mhc import MHC
from neoag_dt.cell.peptide import Peptide
from neoag_dt.cell.allele_score_cache import AlleleScoreCache

_DEFAULT_NAME = "PeptideMHCBindingSimulator"


class pMHC(NamedTuple):
    """ A peptide-MHC complex.
    """
    peptide: Peptide
    hla: MHC


class PeptideMHCBindingSimulator(object):
    """ A class to simulate peptide:MHC binding competition.

    Attributes
    ----------
    binding_scores : neoag_dt.AlleleScoreCache
        The score cache for peptide:MHC binding likelihood (or some other
        binding score such that higher scores are reflective of a higher
        chance to bind)

    name : str
        The name of this object
    """

    def __init__(
            self,
            binding_scores: AlleleScoreCache,
            name: str = _DEFAULT_NAME) -> "PeptideMHCBindingSimulator":

        self.binding_scores = binding_scores
        self.name = name

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def _get_bound_peptides(
            self,
            peptides: Sequence[Peptide],
            hla_allele: MHC,
            num_hla_molecules: int,
            eps: float = 1e-7) -> List[pMHC]:
        """ Select the peptides bound to this HLA """

        if peptides:
            # for each peptide
            #   find its likelihood of binding to this allele
            #
            # these weights will be used for sampling
            peptide_weights = np.array([
                self.binding_scores.get_score(hla_allele, peptide)
                for peptide in peptides
            ])

            # normalize the sampling distribution to sum to 1
            peptide_weights += eps  # avoid peptides with p = 0
            total_weight = np.sum(peptide_weights)
            peptide_weights = peptide_weights / total_weight

            # account for float precision inaccuracy
            if np.any(peptide_weights):
                peptide_weights[0] = (peptide_weights[0] + (1.0 - sum(peptide_weights))).astype(np.float64)

            # for each allele
            #   sample the peptide to which it binds based on the weighted
            if len(peptides) < num_hla_molecules:
                num_hla_molecules = len(peptides)

            bound_peptides = np.random.choice(
                peptides,
                size=num_hla_molecules,
                replace=False,
                p=peptide_weights
            )

            bound_peptides = [
                (peptide, hla_allele) for peptide in bound_peptides
            ]
        else:
            bound_peptides = list()

        return bound_peptides


    def get_bound_peptide_hla_complexes(
            self,
            peptides: Sequence[Peptide],
            hla_counts: Mapping[MHC, int],
            progress_bar: bool = True) -> List[pMHC]:
        """ Simulate a binding competition among `peptides` for the HLA molecules

        Currently, this method models peptide competition among HLA molecules
        of a single type, but not among different HLA molecules. Namely, the
        same number of "starting" peptides are available for each HLA allele.

        More concretely, if a peptide is a "strong, lowly-expressed binder" for
        two HLA alleles, there is not competition for that peptide *between*
        alleles.

        Most concretely, the sampling is *without* replacement for each allele,
        but all peptides are replaced when sampling from the next allele.

        Parameters
        ----------
        peptides : typing.Sequence[neoag_dt.Peptide]
            The list of peptides among which to choose.

        hla_counts : typing.Mapping[neoag_dt.MHC, int]
            The HLAs associated with the current sample, as well as the number
            of molecules of each HLA.

        progress_bar : bool
            Whether to show a progress bar (based on the alleles)

        Returns
        -------
        bound_peptide_hla_complexes : typing.List[neoag_dt.pHMC]
            The peptide:HLA complexes
        """

        it = hla_counts.items()
        if progress_bar:
            it = tqdm.tqdm(it)

        bound_peptides = collection_utils.flatten_lists([
            self._get_bound_peptides(peptides, hla, num_hla_molecules)
            for (hla, num_hla_molecules) in it
        ])

        return bound_peptides
