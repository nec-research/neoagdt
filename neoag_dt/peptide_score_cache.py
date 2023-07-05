""" This module contains a class for storing arbitrary scores associated with
a peptide. In contrast to :obj:`neoag_dt.allele_score_cache.AlleleScoreCache`,
this cache does not store different scores for different HLA alleles.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Mapping

import pandas as pd

from neoag_dt.cell.peptide import Peptide


_DEFAULT_NAME = "PeptideScoreCache"


class PeptideScoreCache(object):
    """ A score cache for peptide sequences

    Attributes
    ----------
    score_cache : typing.Mapping[str, float]
        A map from peptide sequences to a score.

        **N.B.** The keys should be the peptide sequences as strings, not the
        peptide domain objects.

    name : str
        A name for this object
    """
    def __init__(self,
            score_cache:Mapping[str,float],
            name:str=_DEFAULT_NAME) -> "PeptideScoreCache":

        self.score_cache = score_cache
        self.name = name

    def log(self, msg:str, level:int=logging.INFO) -> None:
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def get_score_strings(self, peptide_sequence:str) -> float:
        return self.score_cache.get(peptide_sequence, 0)

    def has_score_strings(self, peptide_sequence:str) -> bool:
        hs = peptide_sequence in self.score_cache
        return hs
        
    def get_score(self, peptide:Peptide) -> float:
        return self.score_cache.get(peptide.sequence, 0)

    def has_score(self, peptide:Peptide) -> bool:
        hs = peptide.sequence in self.score_cache
        return hs
    
    @classmethod
    def construct(klass,
            df_scores:pd.DataFrame,
            peptide_column:str='peptide',
            score_column:str='score',
            name:str=_DEFAULT_NAME) -> "PeptideScoreCache":
        """ Construct a cache using the appropriate columns in `df_scores`

        Parameters
        ----------
        df_scores : pandas.DataFrame
            The data frame

        {peptide,score}_column : str
            The names of the columns containing the respective values

        name : str
            A name for the cache

        Returns
        -------
        peptide_score_cache : neaog_dt.PeptideScoreCache
            The score cache
        """

        keys = df_scores[peptide_column].values
        values = df_scores[score_column].values
        score_cache = dict(zip(keys, values))

        peptide_score_cache = klass(score_cache, name)        
        return peptide_score_cache