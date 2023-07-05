""" This module contains high-level tests regarding the cleavage simulation.
"""
import pytest

import logging
logger = logging.getLogger(__name__)

import pyllars.logging_utils as logging_utils
logging_utils.set_logging_values(logging_level='DEBUG')

from typing import Sequence

import pandas as pd

from neoag_dt import (
    MHC,
    Peptide,
    Protein,
    ProteinCleaver,
    SomaticVariant
)

import neoag_dt.dt_utils as dt_utils
import neoag_dt.test_utils as test_utils

@pytest.fixture
def somatic_variant_name() -> str:
    return '1_1629672_C/F'

@pytest.fixture
def hla_alleles() -> Sequence[MHC]:
    return test_utils.get_hla_alleles()

@pytest.fixture
def df_proteins() -> pd.DataFrame:
    return test_utils.get_proteins_df()

@pytest.fixture
def df_variants() -> pd.DataFrame:
    return test_utils.get_variants_df()

@pytest.fixture
def df_peptide_sequences() -> pd.DataFrame:
    return test_utils.get_peptide_sequences_df()

@pytest.fixture
def somatic_variant() -> SomaticVariant:
    return test_utils.get_somatic_variant(test_utils.get_protein())

@pytest.fixture
def somatic_variant_name() -> str:
    return test_utils.get_somatic_variant(test_utils.get_protein()).name

def test_cleavage(
        somatic_variant: SomaticVariant,
        somatic_variant_name:str,
        hla_alleles:Sequence[MHC],
        df_proteins:pd.DataFrame,
        df_variants:pd.DataFrame,
        df_peptide_sequences:pd.DataFrame) -> None:
    
    # create the maps for peptide extractor
    protein_map = Protein.create_proteins(df_proteins)

    variant_map = SomaticVariant.create_somatic_variants(
        df_variants, protein_map
    )

    peptide_map = Peptide.create_peptides(
        df_peptide_sequences, variant_map
    )
    
    # create random scores for all peptides and alleles
    peptides = list(set(peptide_map.values()))
    presentation_scores = test_utils.create_random_presentation_scores(peptides, hla_alleles)

    # and now the protein cleaver
    protein_cleaver = ProteinCleaver(presentation_scores)
    
    # sample the number of protein molecules based on the expression
    somatic_variant = variant_map[somatic_variant_name]

    protein_mean = somatic_variant.protein.expression_mean
    protein_var = somatic_variant.protein.expression_var
    num_proteins = dt_utils.sample_gamma_poisson(protein_mean, protein_var)

    # and then use the cleaver to select the actual peptides


    selected_peptides = protein_cleaver.cleave(
        peptides, hla_alleles, num_proteins
    )

    assert len(selected_peptides) == num_proteins

def main():
    hla_alleles = test_utils.get_hla_alleles()

    svn = somatic_variant_name()
    
    # load the variants and associated peptide sequences
    df_variants = test_utils.get_variants_df()
    df_peptide_sequences = test_utils.get_peptide_sequences_df()
    df_proteins = test_utils.get_proteins_df()

    test_cleavage(
        somatic_variant,
        somatic_variant_name,
        hla_alleles,
        df_proteins,
        df_variants,
        df_peptide_sequences
    )

if __name__ == '__main__':
    main()