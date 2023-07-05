""" This module contains a class for representing a peptide.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Mapping, Optional, Sequence

import pandas as pd
import pyllars.pandas_utils as pd_utils


class Peptide(object):
    """ A peptide.

    Attributes
    ----------
    somatic_variant : neoag_dt.SomaticVariant
        The variant from which this peptide derives. In case there is not a
        variant associated with this peptide, then this can be null.

    sequence : str
        The peptide sequence
    """
    def __init__(
            self,
            sequence: str,
            somatic_variant: Optional["SomaticVariant"] = None) -> "Peptide":

        self.sequence = sequence
        self.somatic_variant = somatic_variant

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def __hash__(self):
        return hash(self.sequence)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.sequence

    @classmethod
    def create_peptides(
            klass,
            df_peptides: pd.DataFrame,
            variant_map: Mapping[str, "SomaticVariant"],
            sequence_column: str='Mut_peptide',
            variant_column: str = 'Mutation_ID',
            progress_bar: bool = True) -> Mapping[str, "Peptide"]:
        """ Build a mapping from the sequence to the Peptide object

        Parameters
        ----------
        df_peptides : pandas.DataFrame
            A data frame containing all peptides

        variant_map : typing.Mapping[str, neoag_dt.SomaticVariant]
            A mapping from the string name for a variant to the respective
            object

        sequence_column : str
            The name of the column containing the peptide sequence

        progress_bar : bool
            Whether to show a progress bar when constructing the peptides

        Returns
        -------
        peptide_map : typing.Mapping[str, neoag_dt.Peptide]
            A map from each peptide sequence to the respective object

        Notes
        -----
        This function updates all variant objects to know about their
        associated peptide sequences. That is, it has side effects.
        """

        def __from_row(row):
            variant = variant_map.get(row[variant_column])
            p = Peptide(
                sequence=row[sequence_column],
                somatic_variant=variant
            )

            variant.add_peptide(p)
            return p

        peptides = pd_utils.apply(
            df_peptides, __from_row, progress_bar=progress_bar
        )

        peptides = {
            p.sequence : p for p in peptides
        }

        return peptides
