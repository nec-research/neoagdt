""" This module contains a class for representing an HLA allele.
"""
import logging
logger = logging.getLogger(__name__)

import pandas as pd
import lifesci.mhcnames_utils as mhcnames_utils
import pyllars.pandas_utils as pd_utils

from neoag_dt.cell.protein import Protein
from typing import List, Mapping


class MHC(object):
    """ An MHC molecule.

    Attributes
    ----------
    protein : Protein
        The protein object (MHC) associated with this HLA

    name : str
        The name for this HLA.
    """
    def __init__(self,
            protein: Protein,
            name: str) -> "MHC":

        self.name = name
        self.protein = protein

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.name

    @classmethod
    def create_hlas(
            klass,
            df_hlas: pd.DataFrame,
            protein_map: Mapping[str, Protein],
            allele_name_column: str = 'allele_name',
            gene_name_column: str = 'gene_id') -> List["MHC"]:
        """ Create a list of HLAs

        Parameters
        ----------
        df_hlas : pandas.DataFrame
            A data frame containing all HLA information

        protein_map : typing.Mapping[str, neoag_dt.Protein]
            A map from each protein name to the respective object

        {allele,protein}_name_column : str
            The names of the respective columns in the data frame
        """

        def __from_row(row):
            hla_name = row[allele_name_column]
            gene_name = row[gene_name_column]
            protein = protein_map[gene_name]
            h = MHC(protein=protein, name=hla_name)
            return h

        hlas = pd_utils.apply(df_hlas, __from_row)
        return hlas
