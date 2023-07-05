""" This module contains a class for simulating protein cleavage for selecting
peptides.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List, Sequence
import numpy as np

from neoag_dt.cell.mhc import MHC
from neoag_dt.cell.peptide import Peptide
from neoag_dt.cell.allele_score_cache import AlleleScoreCache


_DEFAULT_NAME = "ProteinCleaver"


class ProteinCleaver(object):
    """ A simulation of protein cleavage.

    Attributes
    ----------
    cleavage_scores : neoag_dt.AlleleScoreCache
        The score cache for peptide:MHC presentation likelihood

    name : str
        The name of this object
    """
    def __init__(
            self,
            cleavage_scores: AlleleScoreCache,
            name: str = _DEFAULT_NAME) -> "ProteinCleaver":

        self.name = name
        self.cleavage_scores = cleavage_scores

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def cleave(
            self,
            peptides: Sequence[Peptide],
            hla_alleles: Sequence[MHC],
            num_proteins: int = 1) -> List[Peptide]:
        """ Simulate protein cleavage for selecting among `peptides`

        The cleavage simulation uses the cleavage prediction for each
        peptide across all alleles to create a weighted distribution. This
        distribution is then sampled from for `num_proteins` times to determine
        which peptide was created during each cleavage step.

        Parameters
        ----------
        peptides : typing.Sequence[neoag_dt.Peptide]
            The list of peptides among which to choose. Presumably, this list
            is based on a sliding window around a variant, or for all positions
            downstream of an indel, or similar.

        hla_alleles : typing.Sequence[neoag_dt.MHC]
            The list of HLAs associated with the current sample.

        num_proteins : int
            The number of times to simulate the cleavage operation

        Returns
        -------
        selected_peptides : typing.List[neoag_dt.Peptide]
            The peptides selected by each cleavage simulation
        """

        # these are the weights for sampling
        peptide_weights = np.array([
            self.cleavage_scores.get_maximum_peptide_score(pep, hla_alleles)
                for pep in peptides
        ])

        # normalize the sampling distribution to sum to 1
        total_weight = np.sum(peptide_weights)
        peptide_weights = peptide_weights / total_weight

        # for each protein molecule
        #   sample the peptide based on the weighted 
        selected_peptides = np.random.choice(
            peptides,
            size=num_proteins,
            replace=True,
            p=peptide_weights
        )

        return selected_peptides
