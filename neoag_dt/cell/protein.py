""" This module contains a class for representing a protein.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Mapping, Optional

import pandas as pd
import pyllars.pandas_utils as pd_utils

from neoag_dt.cell.somatic_variant import SomaticVariant


class Protein(object):
    """ A protein.

    These objects are called "proteins" for semantic reasons. In practice,
    they are likely to be RNA transcripts, or even just genes, since
    RNA-seq is the most common way to actually measure gene expression.

    Attributes
    ----------
    expression_{mean,var} : float
        The mean and variance of the expression (e.g., TPM) for this "protein".

        If possible, it is better if this is the value for the actual
        variant to which this "protein" will eventually be linked.

        **N.B.** This value should not be log-transformed for use in these
        simulations.

    name : str
        The name for this protein.

    variants : List[neoag_dt.SomaticVariant]
        The variants associated with this protein
    """
    def __init__(
            self,
            name: str,
            expression_mean: Optional[float] = None,
            expression_var: Optional[float] = None) -> "Protein":

        self.name = name
        self.expression_mean = expression_mean
        self.expression_var = expression_var
        self.variants = []

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def add_variant(self, variant:SomaticVariant) -> "Protein":
        self.variants.append(variant)
        return self

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.name

    @classmethod
    def create_proteins(
            klass,
            df_proteins: pd.DataFrame,
            name_column: str = 'gene_id',
            expression_mean_column: Optional[str] = 'FPKM',
            expression_var_column: Optional[str] = 'FPKM_VAR',
            default_mean_value: Optional[float] = 0,
            default_var_value: Optional[float] = 1,
            progress_bar: bool = True) -> Mapping[str, "Protein"]:
        """ Build a mapping from the protein name to the Protein object

        **N.B.** This function extracts the unique protein names and only
        includes one entry per protein.
        
        Parameters
        ----------
        df_proteins : pandas.DataFrame
            A data frame containing all proteins

        name_column : str
            The name of the column containing the protein name

        expression_mean_column : typing.Optional[str]
            The columns containing the expression (such as TPM) mean

        expression_var_column : typing.Optional[str]
            The columns containing the expression (such as TPM) variance

        default_mean_value : typing.Optional[float]
            Default values for the respective mean values

        default_var_value : typing.Optional[float]
            Default values for the respective variance values

        progress_bar : bool
            Whether to show a progress bar when constructing the peptides

        Returns
        -------
        protein_map : typing.Mapping[str, neoag_dt.Protein]
            A map from each protein name to the respective object
        """

        def __from_row(row):
            v = Protein(
                name=row[name_column],
                expression_mean=row.get(expression_mean_column, default_mean_value),
                expression_var=row.get(expression_var_column, default_var_value)
            )
            return v

        df_proteins = df_proteins.drop_duplicates(subset=[name_column])

        if expression_var_column is not None:
            msg = ("[Protein.create_proteins] Replacing all non-positive "
                "expression variances with {:.2f}\n\nN.B. This makes a copy of "
                "the data frame.".format(default_var_value))

            df_proteins = df_proteins.copy()
            m_non_positive = (df_proteins[expression_var_column] <= 0)
            df_proteins.loc[m_non_positive, expression_var_column] = default_var_value

        proteins = pd_utils.apply(
            df_proteins, __from_row, progress_bar=progress_bar
        )

        proteins = {
            p.name: p for p in proteins
        }

        return proteins
