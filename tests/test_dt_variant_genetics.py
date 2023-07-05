""" This module contains high-level tests for the "genetics" aspect of the
digital twin for variants. Specifically, "genetics" means the DNA, RNA, protein
molecule steps in the process.
"""
import pytest

import logging
logger = logging.getLogger(__name__)

import pyllars.logging_utils as logging_utils
logging_utils.set_logging_values(logging_level='DEBUG')

import pandas as pd

from neoag_dt import (
    Protein, SomaticVariant, GeneticSimulator
)
import neoag_dt.test_utils as test_utils


@pytest.fixture
def somatic_variant(protein:Protein) -> SomaticVariant:
    return test_utils.get_somatic_variant(protein)

@pytest.fixture
def protein() -> Protein:
    return test_utils.get_protein()

@pytest.fixture
def df_proteins() -> pd.DataFrame:
    return test_utils.get_proteins_df()

@pytest.fixture
def df_variants() -> pd.DataFrame:
    return test_utils.get_variants_df()

def test_genetic_simulator(somatic_variant:SomaticVariant) -> None:
    simulator = GeneticSimulator()
    result = simulator.simulate(somatic_variant)

    assert isinstance(result.variant_in_dna, bool)
    assert isinstance(result.protein_count, int)

def test_genetics(
        df_proteins:pd.DataFrame,
        df_variants:pd.DataFrame) -> None:
    protein_map = Protein.create_proteins(df_proteins)

    variant_map = SomaticVariant.create_somatic_variants(
        df_variants, protein_map
    )

    # for each variant, we need to check RNA and protein count
    simulator = GeneticSimulator()
    results = simulator.simulate_all(variant_map)

    assert len(results) == len(variant_map)

def main():
    
    # load some simple domain objects for initial tests
    protein = test_utils.get_protein()
    somatic_variant = test_utils.get_somatic_variant(protein)

    # and actual data for later tests
    df_variants = test_utils.get_variants_df()
    

    test_genetic_simulator(somatic_variant)
    test_genetics(df_variants)

if __name__ == '__main__':
    main()