""" This module tests the entire process of simulating a cancer cell.
"""
import pytest

import logging
logger = logging.getLogger(__name__)

import pyllars.logging_utils as logging_utils
logging_utils.set_logging_values(logging_level='DEBUG')

from typing import Sequence

import numpy as np
import pandas as pd

from neoag_dt import (
    CellFactory, GeneticSimulator, MHC, Peptide,
    PeptideMHCBindingSimulator, Protein, ProteinCleaver,
    SomaticVariant
)
import neoag_dt.test_utils as test_utils

@pytest.fixture
def hla_alleles() -> Sequence[MHC]:
    return test_utils.get_hla_alleles()

@pytest.fixture
def df_variants() -> pd.DataFrame:
    return test_utils.get_variants_df()

@pytest.fixture
def df_peptide_sequences() -> pd.DataFrame:
    return test_utils.get_peptide_sequences_df()

@pytest.fixture
def df_proteins() -> pd.DataFrame:
    return test_utils.get_proteins_df()

@pytest.fixture
def df_hlas() -> pd.DataFrame:
    return test_utils.get_hlas_df()

def test_hla_loading(
        df_proteins:pd.DataFrame,
        df_hlas:pd.DataFrame) -> None:

    # create the domain objects
    protein_map = Protein.create_proteins(df_proteins)
    hlas = MHC.create_hlas(df_hlas, protein_map)

    assert len(hlas) == 6

def test_cell_factory(
        df_proteins:pd.DataFrame,
        df_variants:pd.DataFrame,
        df_peptide_sequences:pd.DataFrame,
        hla_alleles:Sequence[MHC]) -> None:

    # create the domain objects
    protein_map = Protein.create_proteins(df_proteins)

    variant_map = SomaticVariant.create_somatic_variants(
        df_variants, protein_map
    )

    peptide_map = Peptide.create_peptides(
        df_peptide_sequences, variant_map
    )

    # create random scores for all peptides and alleles
    peptides = peptide_map.values()
    presentation_scores = test_utils.create_random_presentation_scores(peptides, hla_alleles)
    binding_scores = test_utils.create_random_binding_scores(peptides, hla_alleles)

    # and the object to perform the cleavage
    protein_cleaver = ProteinCleaver(presentation_scores)

    # for each variant, we need to check DNA, RNA, and protein count
    genetic_simulator = GeneticSimulator()

    # and binding competition
    binding_simulator = PeptideMHCBindingSimulator(binding_scores)


    ###
    # Finally, create and test the cell factory
    ###
    cell_factory = CellFactory(
        genetic_simulator=genetic_simulator,
        protein_cleaver=protein_cleaver,
        binding_simulator=binding_simulator,
        presentation_scores=presentation_scores
    )

    cell = cell_factory.create_cell(
        variant_map, hla_alleles, name='cell'
    )
    
    expected_num_peptides = np.sum([
        r.protein_count for r in cell.genetic_simulation_results.values()
    ])

    actual_num_peptides = len(cell.selected_peptides)

    assert actual_num_peptides == expected_num_peptides

    num_hlas = np.sum([
        cell.hla_counts[hla] for hla in hla_alleles
    ])

    assert num_hlas > 0

    # There can be less peptides than MHC molecules to which they could bind
    assert len(cell.bound_peptide_hlas) <= num_hlas

    # make sure that we have a feasible number of peptides presented

    # In principle, it could happen that all of the peptides are presented,
    # but that is exceedingly unlikely for more than a few peptides. Assume
    # that even if this happens, it is still something to look at.
    #
    # Thus, use "<" rather than "<="
    assert len(cell.presented_peptide_hlas) < num_hlas

    all_pmhc = [
        {
            'peptide': pmhc.peptide.sequence,
            'hla': pmhc.hla.name,
            'cell': cell.name
         } for pmhc in cell.presented_peptide_hlas
    ]

    df_pmhc = pd.DataFrame(all_pmhc)
    return df_pmhc

def test_cell_factory_empty_cell(
        df_proteins: pd.DataFrame,
        hla_alleles: Sequence[MHC]) -> None:
    # create the domain objects and call functions with empty peptide and variant dataframes
    protein_map = Protein.create_proteins(df_proteins)

    variant_map = SomaticVariant.create_somatic_variants(
        pd.DataFrame(), protein_map
    )

    peptide_map = Peptide.create_peptides(
        pd.DataFrame(), variant_map
    )
    peptides = peptide_map.values()
    presentation_scores = test_utils.create_random_presentation_scores(peptides, hla_alleles)
    binding_scores = test_utils.create_random_binding_scores(peptides, hla_alleles)

    # There should be no scores
    assert len(binding_scores) == 0
    assert len(presentation_scores) == 0

    protein_cleaver = ProteinCleaver(presentation_scores)
    genetic_simulator = GeneticSimulator()
    binding_simulator = PeptideMHCBindingSimulator(binding_scores)

    cell_factory = CellFactory(
        genetic_simulator=genetic_simulator,
        protein_cleaver=protein_cleaver,
        binding_simulator=binding_simulator,
        presentation_scores=presentation_scores
    )

    cell = cell_factory.create_cell(
        variant_map, hla_alleles, name='cell'
    )

    # There shouldn't be any peptides selected
    assert len(cell.selected_peptides) == 0

    # There shouldn't be any presented pMHCs
    assert len(cell.bound_peptide_hlas) == 0

def main():
    hla_alleles = test_utils.get_hla_alleles()
    
    # load the variants and associated peptide sequences
    df_variants = test_utils.get_variants_df()
    df_peptide_sequences = test_utils.get_peptide_sequences_df()
    df_proteins = test_utils.get_proteins_df()
    df_hlas = test_utils.get_hlas_df()

    test_hla_loading(df_proteins, df_hlas)

    test_cell_factory(df_proteins, df_variants, df_peptide_sequences, hla_alleles)

    # test for ccells without peptides
    test_cell_factory_empty_cell(df_proteins, hla_alleles)

if __name__ == '__main__':
    main()