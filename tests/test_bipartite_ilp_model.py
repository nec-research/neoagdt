"""This module tests the BipartiteVaccineDesignModel class.
"""
import pytest

import logging
logger = logging.getLogger(__name__)

import pyllars.logging_utils as logging_utils
logging_utils.set_logging_values(logging_level='DEBUG')

import mip

from neoag_dt import BipartiteVaccineDesignModel
from neoag_dt.optimization import optimization_utils
import neoag_dt.test_utils as test_utils
from neoag_dt.cli.optimize_vaccine_ilp import get_weights


def test_approximate_response_scores():
    ret = test_utils.get_ilp_model_input()
    optim_config, peptides, df_peptides, df_cells, num_cells, p_response_factor = ret

    df_scores = optimization_utils.approximate_response_scores(
        df_cells,
        peptides,
        num_cells,
        optim_config
    )

    dfs_map = test_utils.get_distance_from_self_map()

    for cell in df_cells['cell_ids'].unique():
        df_scores_cell = df_scores[df_scores['cell_id'] == cell]
        presented_peptides = df_cells[df_cells['cell_ids'] == cell]['presented_peptides'].tolist()
        for peptide in peptides:
            p_response = df_scores_cell[df_scores_cell['vaccine_element'] == peptide]['p_response']
            if peptide in presented_peptides:
                presentation_count = presented_peptides.count(peptide)
                expected_p_value = min(1, presentation_count * p_response_factor) - 1e-7
                expected_p_value *= dfs_map[peptide]
                assert p_response.iloc[0] == expected_p_value
            else:
                assert p_response.iloc[0] == 1e-7


def test_get_weights():
    ret = test_utils.get_ilp_model_input()
    optim_config, peptides, df_peptides, df_cells, num_cells, p_response_factor = ret

    peptides_weights_map = get_weights(
        df_peptides, 'Mut_peptide', None
    )

    for peptide in peptides_weights_map.keys():
        assert peptides_weights_map[peptide] == 1


def test_bipartite_ilp_model():
    ret = test_utils.get_ilp_model_input(criterion='MinSum')
    optim_config, peptides, df_peptides, df_cells, num_cells, p_response_factor = ret

    df_scores = optimization_utils.approximate_response_scores(
        df_cells,
        peptides,
        num_cells,
        optim_config
    )

    peptides_weights_map = get_weights(
        df_peptides, 'Mut_peptide', None
    )

    ilp_model = BipartiteVaccineDesignModel(
        config=optim_config,
        peptides=peptides,
        df_cells=df_cells,
        df_scores=df_scores,
        peptide_weights_map=peptides_weights_map,
        vaccine_elements="peptides",
        seed=42
    )

    for peptide in ilp_model.peptide_index_map.keys():
        assert peptide in peptides

    assert ilp_model.num_peptides == len(peptides)

    for cell in ilp_model.cell_index_map.keys():
        assert cell in df_cells['cell_ids'].unique()

    assert ilp_model.num_cells == len(df_cells['cell_ids'].unique())

    ilp_model.build_model()

    assert len(ilp_model.x_peptide) == len(peptides)
    for x in ilp_model.x_peptide:
        assert x.var_type == mip.BINARY

    assert len(ilp_model.x_cell) == len(df_cells['cell_ids'].unique())
    for x in ilp_model.x_cell:
        assert x.var_type == mip.CONTINUOUS

    ilp_model.optimize()

    df_selected_peptides = ilp_model.get_selected_peptides()

    for rep in df_cells['repetition'].unique():
        assert rep in df_selected_peptides['repetition'].unique()
        selected_peptides_for_rep = len(
            df_selected_peptides[df_selected_peptides['repetition'] == rep]
        )
        assert selected_peptides_for_rep <= optim_config['budget']
