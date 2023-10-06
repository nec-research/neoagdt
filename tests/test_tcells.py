""" This module contains tests for generating the tumor microenvironment.
"""
import pytest

import logging
logger = logging.getLogger(__name__)

import pyllars.logging_utils as logging_utils
logging_utils.set_logging_values(logging_level='DEBUG')

from typing import List, Sequence

from neoag_dt import (
    TCellActivationSimulator, TCellResponseSimulator, Peptide, PeptideVaccine, TCell, pMHC, Cell
)
import neoag_dt.test_utils as test_utils

@pytest.fixture
def peptide_vaccine() -> PeptideVaccine:
    return test_utils.get_pan_pool_1_vaccine()

@pytest.fixture
def peptides() -> List[Peptide]:
    return test_utils.get_peptides()

@pytest.fixture
def activated_t_cells(peptides:Sequence[Peptide]) -> List[TCell]:
    return test_utils.get_activated_t_cells(peptides)

@pytest.fixture
def tumor_cells(peptides:Sequence[Peptide]) -> List[Cell]:
    hla = test_utils.get_hla_allele()
    return test_utils.get_tumor_cells(peptides, hla)

def test_t_cell_activation(peptide_vaccine:PeptideVaccine) -> None:

    dfs_cache = test_utils.create_random_dfs_scores(peptide_vaccine.peptides)

    t_cell_simulator = TCellActivationSimulator(dfs_cache=dfs_cache)
    activated_t_cells = t_cell_simulator.simulate(peptide_vaccine)

    # it is very unlikely that we will have no activated T cells
    assert len(activated_t_cells) > 0

def test_t_cell_responses(
        activated_t_cells:Sequence[TCell],
        tumor_cells:Sequence[Cell]) -> None:

    peptides = {p
        for c in tumor_cells
            for p in c.presented_peptides
    }
    response_cache = test_utils.create_random_response_scores(peptides)

    response_simulator = TCellResponseSimulator(activated_t_cells, response_cache)

    cells_with_response = response_simulator.simulate(tumor_cells)

    assert len(cells_with_response) > 0



def main():
    # given vaccine elements
    # right now, these are just peptides
    vaccine = test_utils.get_pan_pool_1_vaccine()

    # for each vaccine element
    # does this activate T cells?
    # if so, how many?
    test_t_cell_activation(vaccine)


    peptides = test_utils.get_peptides()
    activated_t_cells = test_utils.get_activated_t_cells(peptides)
    hla = test_utils.get_hla_allele()
    tumor_cells = test_utils.get_tumor_cells(peptides, hla)

    test_t_cell_responses(activated_t_cells, tumor_cells)

if __name__ == '__main__':
    main()