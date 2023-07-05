""" This module contains a class for storing arbitrary scores associated with
(HLA, peptide) pairs. In contrast to :obj:`neoag_dt.peptide_score_cache.PeptideScoreCache`,
this cache requires an HLA to be associated with each score.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List, Mapping, Sequence, Tuple

import pandas as pd

allele_ve_type = Tuple[str,str]
score_cache_type = Mapping[allele_ve_type, float]

from neoag_dt.cell.mhc import MHC
from neoag_dt.cell.peptide import Peptide

_DEFAULT_NAME = "AlleleScoreCache"


class AlleleScoreCache:
    """ A score cache for (HLA, peptide) pairs.

    Attributes
    ----------
    score_cache : typing.Mapping[typing.Tuple[str,str], float]
        A map from (HLA allele name, peptide sequence) pairs to a score.

        **N.B.** The order of the scores must be allele, then peptide sequence.
        Additionaly, note that the keys should be the allele name and peptide
        sequences as strings, not the allele and peptide domain objects.

    name : str
        A name for this object
    """

    def __init__(
            self,
            score_cache: score_cache_type,
            name: str = _DEFAULT_NAME) -> "AlleleScoreCache":

        self.score_cache = score_cache
        self.name = name

    def log(self, msg: str, level: int = logging.INFO) -> None:
        """ Log `msg` using `level` using the module-level logger """
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def get_score_strings(self, allele_name: str, peptide_sequence: str) -> float:
        """ Get scores given allele and peptide"""
        return self.score_cache.get((allele_name, peptide_sequence), 0)

    def has_score_strings(self, allele_name: str, peptide_sequence: str) -> bool:
        """Check cache"""
        has = (allele_name, peptide_sequence) in self.score_cache
        return has

    def get_score(self, allele: MHC, peptide: Peptide) -> float:
        """Get score given MHC allele and peptide"""
        return self.score_cache.get((allele.name, peptide.sequence), 0)

    def has_score(self, allele: MHC, peptide: Peptide) -> bool:
        """Check score is in cache"""
        has = (allele.name, peptide.sequence) in self.score_cache
        return has

    def get_minimum_score(self) -> float:
        """Get min score"""
        minimum_score = min(self.score_cache.values())
        return minimum_score

    def get_maximum_score(self) -> float:
        """Get max score"""
        maximum_score = max(self.score_cache.values())
        return maximum_score

    def get_maximum_peptide_score(self,
                                  peptide: Peptide,
                                  alleles: Sequence[MHC]) -> float:
        """ Retrieve the maximum score for `peptide_sequence` among all of the
        alleles in `allele_names`

        Parameters
        ----------
        peptide : Peptide
            The peptide

        alleles : typing.Sequence[MHC]
            The alleles to consider

        Returns
        -------
        maximum_score : float
            The maximum score for the sequence among all of the given alleles
        """
        scores = [
            self.get_score(allele, peptide)
                for allele in alleles
        ]

        max_score = max(scores)
        return max_score

    def get_all_allele_names(self) -> List[str]:
        """Get all ellele names"""
        all_allele_names = {
            k[0] for k in self.score_cache.keys()
        }
        return all_allele_names

    def __len__(self):
        return len(self.score_cache)

    @classmethod
    def construct(
            klass,
            df_scores: pd.DataFrame,
            allele_column: str = 'allele',
            peptide_column: str = 'peptide',
            score_column: str = 'score',
            name: str = _DEFAULT_NAME) -> "AlleleScoreCache":
        """ Construct a cache using the appropriate columns in `df_scores`

        Parameters
        ----------
        df_scores : pandas.DataFrame
            The data frame

        {allele,peptide,score}_column : str
            The names of the columns containing the respective values

        name : str
            A name for the cache

        Returns
        -------
        allele_score_cache : neaog_dt.AlleleScoreCache
            The score cache
        """

        alleles = df_scores[allele_column]
        vaccine_elements = df_scores[peptide_column]

        keys = list(zip(alleles, vaccine_elements))
        values = df_scores[score_column]
        score_cache = dict(zip(keys, values))

        allele_score_cache = klass(score_cache, name)
        return allele_score_cache
