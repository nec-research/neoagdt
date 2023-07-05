""" This module contains a class for representing a somatic variant.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Mapping, Optional, Sequence, Union

import pandas as pd
import pyllars.pandas_utils as pd_utils

from neoag_dt.cell.peptide import Peptide


class SomaticVariant(object):
    """ A somatic variant.

    Attributes
    ----------
    dna_ref_count : int
        The count of DNA (presumably, WES) sequencing reads which contain the
        reference nucleotide for this variant.

    dna_alt_count : int
        The count of DNA (presumably, WES) sequencing reads which contain the
        variant ("alternative") nucleotide for this variant.

    rna_ref_count : int
        The count of RNA sequencing reads which contain the reference
        nucleotide for this variant.

    rna_alt_count : int
        The count of RNA sequencing reads which contain the variant
        ("alternative") nucleotide for this variant.

    vaf : typing.Optional[float]
        The VAF; if missing, it is computed.

    protein : typing.Optional[neoag_dt.Protein]
        The source protein (transcript, gene, etc., depending on the level at
        which gene expression was actually measure in the real experiments)
        for this variant.

    name : str
        A name for this variant

    peptides : typing.List[neoag_dt.Peptide]
        The peptides associated with this variant
    """ 

    def __init__(
            self,
            dna_ref_count: Optional[int],
            dna_alt_count: Optional[int],
            rna_ref_count: Optional[int],
            rna_alt_count: Optional[int],
            vaf: Optional[float] = None,
            protein: Optional["Protein"] = None,
            name: Optional[str] = None) -> "SomaticVariant":

        self.name = name
        self.dna_ref_count = dna_ref_count
        self.dna_alt_count = dna_alt_count
        self.rna_ref_count = rna_ref_count
        self.rna_alt_count = rna_alt_count
        self.vaf = vaf
        self.protein = protein
        self.peptides = []

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def add_peptide(self, peptide: Peptide) -> "SomaticVariant":
        self.peptides.append(peptide)
        return self

    def calculate_vaf_dna(self):
        if self.vaf is not None:
            return self.vaf
        else:
            dna_count = self.dna_ref_count + self.dna_alt_count
            if dna_count == 0:
                af = 0.0
            else:
                af = self.dna_alt_count / dna_count
            return af

    def calculate_vaf_rna(self):
        if self.vaf is not None:
            return self.vaf
        else:
            rna_count = self.rna_ref_count + self.rna_alt_count
            if rna_count == 0:
                af = 0.0
            else:
                af = self.rna_alt_count / rna_count
            return af

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not(self == other)

    def __str__(self):
        return self.name

    @classmethod
    def create_somatic_variants(
            klass,
            df_variants: pd.DataFrame,
            protein_map: Mapping[str, "Protein"],
            name_column: str = 'Mutation_ID',
            dna_ref_count_column: Optional[str] = 'WXS_tumor_depth_ref',
            dna_alt_count_column: Optional[str] = 'WXS_tumor_depth_alt',
            rna_ref_count_column: Optional[str] = 'RNA_tumor_depth_ref',
            rna_alt_count_column: Optional[str] = 'RNA_tumor_depth_alt',
            vaf_column: Optional[str] = None,
            gene_column: str = 'Gene_ID',
            progress_bar: bool = True) -> Mapping[str, "SomaticVariant"]:
        """ Build a mapping from the variant name to the SomaticVariant object
        
        Parameters
        ----------
        df_variants : pandas.DataFrame
            The data frame containing the variant information

        protein_map : typing.Mapping[str, neoag_dt.Protein]
            A mapping from protein names to Protein objects

        name_column : str
            The name of the column in `df_variants` which contains the name
            of each variant. These names must match those found in 
            `df_variant_peptide_sequences`.

        dna_ref_count_column : typing.Optional[str]
            The columns containing the respective counts in `df_variants`

        dna_alt_count_column : typing.Optional[str]
            The columns containing the respective counts in `df_variants`

        rna_ref_count_column : typing.Optional[str]
            The columns containing the respective counts in `df_variants`

        rna_alt_count_column : typing.Optional[str]
            The columns containing the respective counts in `df_variants`

        vaf_column : typing.Optional[str]
            The columns containing the respective counts in `df_variants`

        gene_column : str
            The column containing the source protein (transcript, gene, etc.)
            for this variant

        progress_bar : bool
            Whether to show a progress bar when constructing the peptides

        Returns
        -------
        somatic_variant_map : typing.Mapping[str, neoag_dt.SomaticVariant]
            A map from each somatic variant name to the respective object

        Notes
        -----
        This function updates all protein objects to know about their somatic
        variants. That is, it has side effects.
        """

        def __from_row(row):
            protein = protein_map.get(row[gene_column])

            if vaf_column is None:
                # the vaf can be None,
                # in that case it is computed in the `SomaticVariant` class using
                # the `dna_ref_count_column`, `dna_alt_count_column`, `rna_ref_count_column` and
                # `rna_alt_count_column`

                assert (
                    dna_ref_count_column is not None and
                    dna_alt_count_column is not None and
                    rna_ref_count_column is not None and
                    rna_alt_count_column is not None
                ), "If `vaf` is None, `var_dna_ref`, `var_dna_alt`, `var_rna_ref`" \
                   "and `var_rna_alt` cannot be None"

                v = SomaticVariant(
                    dna_ref_count=row[dna_ref_count_column],
                    dna_alt_count=row[dna_alt_count_column],
                    rna_ref_count=row[rna_ref_count_column],
                    rna_alt_count=row[rna_alt_count_column],
                    protein=protein,
                    name=row[name_column],
                    vaf=None
                )
            else:
                v = SomaticVariant(
                    dna_ref_count=None,
                    dna_alt_count=None,
                    rna_ref_count=None,
                    rna_alt_count=None,
                    protein=protein,
                    name=row[name_column],
                    vaf=row[vaf_column]
                )

            protein.add_variant(v)

            return v

        somatic_variants = pd_utils.apply(
            df_variants, __from_row, progress_bar=progress_bar
        )

        somatic_variants = {
            v.name : v for v in somatic_variants
        }

        return somatic_variants
