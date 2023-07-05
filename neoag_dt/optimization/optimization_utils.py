""" This module contains helpers for running the optimization models.
"""
import logging
logger = logging.getLogger(__name__)

import pandas as pd
import math

import pyllars.collection_utils as collection_utils

from typing import Mapping, Sequence, Union


SIMULATION_COLS = [
    'simulation_name',
    'repetition'
]

###
# Determining which simulations were generated
#
#    This functionality is here because it needs to match with
#    `SIMULATION_COLS`, which we later use when optimizing peptide
#    selection.
###
def _get_simulation_group_vals(simulation_info):
    name = simulation_info['simulation_name']
    num_repetitions = simulation_info['num_repetitions']
    num_cells = simulation_info['simulation_num_cells']

    ret = [
        (name, repetition, num_cells) for repetition in range(num_repetitions)
    ]
    
    return ret


def get_simulation_group_vals(simulation_config):
    simulations = simulation_config['simulations']

    sim_info = [
        _get_simulation_group_vals(simulation_info)
            for simulation_info in simulations
    ]

    sim_info = collection_utils.flatten_lists(sim_info)
    return sim_info


def approximate_response_scores(
        df_cells: pd.DataFrame,
        peptides: Sequence[str],
        num_cells: int,
        optimization_config: Mapping[str, Union[str, float, int]],
        epsilon: float = 1e-7) -> pd.DataFrame:
    """ Approximate probabilities of no response logP(R=- | v_i, c_j).

    Assumption:
    If a cancer cell c_j presents a peptide v_i N times, including vaccine element v_i
    in a vaccine guarantees an immune response defined by:

    `P(R=+ | v_i, c_j) = (min(1, N * p_response_factor) - epsilon) * distance_from_self(v_i)`
    """
    data = {
        'cell_id': [],
        'vaccine_element': [],
        'peptide_len': [],
        'p_response': [],
        'log_p_response': [],
        'log_p_no_response': []
    }

    p_response_factor = optimization_config['p_response_factor']
    if p_response_factor is None:
        p_response_factor = get_adaptive_p_response_factor(df_cells)

    dfs_map = get_dfs(optimization_config)

    for cell_id in range(num_cells):
        presented_peptides = df_cells[df_cells['cell_ids'] == cell_id]['presented_peptides'].tolist()
        for peptide in peptides:
            data['cell_id'].append(cell_id)
            data['vaccine_element'].append(peptide)
            data['peptide_len'].append(len(peptide))
            if peptide in presented_peptides:
                presentation_count = presented_peptides.count(peptide)
                p_response = min(1, presentation_count * p_response_factor) - epsilon
                if dfs_map:  # distance-from-self contribution to P(R=+ | v_i, c_j)
                    p_response *= dfs_map[peptide]
                p_no_response = 1 - p_response
                data['p_response'].append(p_response)
                data['log_p_response'].append(math.log(p_response))
                data['log_p_no_response'].append(math.log(p_no_response))
            else:
                data['p_response'].append(epsilon)
                data['log_p_response'].append(math.log(epsilon))
                data['log_p_no_response'].append(math.log(1 - epsilon))

    df_scores = pd.DataFrame(data=data)
    return df_scores


def get_cell_peptide_map(df_cells: pd.DataFrame, peptides: Sequence[str]) -> Mapping:
    """ Maps each cell to the set of presented peptides.
    """
    cell_ids = df_cells['cell_ids'].unique()
    cell_peptide_map = {cell_id: [] for cell_id in cell_ids}

    for cell_id in cell_ids:
        peptides_on_cell = df_cells[df_cells['cell_ids'] == cell_id]['presented_peptides'].unique()
        for peptide in peptides:
            if peptide in peptides_on_cell:
                cell_peptide_map[cell_id].append(peptide)

    return cell_peptide_map


def get_index_and_reverse_map(items):
    index_map = {
        c: i for i, c in enumerate(items)
    }
    reverse_index_map = {
        i: c for c, i in index_map.items()
    }

    return index_map, reverse_index_map


###
# Processing base scores for creating MIP
###
def get_peptide_cell_scores_matrix(df_scores, peptide_index_map, cell_index_map):
    peptide_indices = df_scores['vaccine_element'].map(peptide_index_map)
    allele_indices = df_scores['cell_id'].map(cell_index_map)
    keys = zip(peptide_indices, allele_indices)

    # and values will be the log probability of no response
    values = df_scores['log_p_no_response']
    peptide_cell_scores = dict(zip(keys, values))

    return peptide_cell_scores


def extract_best_peptides(df: pd.DataFrame, budget: int) -> pd.DataFrame:
    """ Inputs the results of multiple cell population simulations
    (multiple repetitions and multiple population sizes) and
    extract the peptides which were selected most times.
    """
    peptides_weights_map = dict(zip(df['peptide'], df['weight']))

    # count how many times the peptides were selected
    df = df.groupby(['peptide']).size().reset_index(name='counts')
    df = df.sort_values(by=['counts'], ascending=False)
    peptides_counts_map = dict(zip(df['peptide'], df['counts']))

    vaccine_elements = []
    total_weight = 0

    # select peptides which were selected most times
    for peptide in df['peptide'].tolist():
        peptide_weight = peptides_weights_map[peptide]
        if total_weight + peptide_weight <= budget:
            vaccine_elements.append(peptide)
            total_weight += peptide_weight

    vaccine_df = pd.DataFrame({'peptide': vaccine_elements})
    vaccine_df['counts'] = vaccine_df['peptide'].apply(lambda x: peptides_counts_map[x])
    vaccine_df['weight'] = vaccine_df['peptide'].apply(lambda x: peptides_weights_map[x])

    return vaccine_df


def get_dfs(optimization_config: Mapping[str, Union[str, float, int]]):
    """ Load distance-from-self.
    If optimization_config['distance_from_self_scores'] is None,
    we return a None dictionary.
    """
    if optimization_config['distance_from_self_scores'] is not None:
        msg = f"loading distance-from-self scores: {optimization_config['distance_from_self_scores']}"
        logger.info(msg)
        df_dfs = pd.read_csv(optimization_config['distance_from_self_scores'])
        dfs_map = dict(zip(
            df_dfs[optimization_config['peptide_sequence_column']],
            df_dfs[optimization_config['distance_from_self_column']]
        ))
    else:
        msg = f"Not using distance-from-self."
        logger.info(msg)
        dfs_map = None

    return dfs_map


def aggregate_by_mutations(
        df_scores: pd.DataFrame,
        df_peptides: pd.DataFrame,
        optimization_config: Mapping[str, Union[str, float, int]]
    ) -> pd.DataFrame:
    """ This function is only called when the optimization is run over
    mutations. It computes the log-probability of no response for a
    given mutation. This is computed by aggregating the peptide-specific:
    log P(R =  - | v_i, c_j) = log \prod_{k=1}^{O^i} P(R =  - | p^i_k, c_j),
    where v_i is a mutation, p^i_k is a minimal epitope, O^i is the number of
    minimal epitopes associated to v_i.
    """
    pep_mutation_map = dict(zip(
        df_peptides[optimization_config['peptide_sequence_column']],
        df_peptides[optimization_config['mutation_id_column']],
    ))

    df_scores['mutation'] = df_scores['vaccine_element'].apply(lambda x: pep_mutation_map[x])
    df_scores = df_scores.groupby(['cell_id', 'mutation'])['log_p_no_response'].sum().reset_index()

    df_scores['p_no_response'] = df_scores['log_p_no_response'].apply(lambda x: math.exp(x))
    df_scores['p_response'] = df_scores['p_no_response'].apply(lambda x: 1 - x)
    df_scores['log_p_response'] = df_scores['p_response'].apply(lambda x: math.log(x))
    df_scores = df_scores.rename(columns={'mutation': 'vaccine_element'})

    return df_scores


def get_adaptive_p_response_factor(df_cells: pd.DataFrame) -> float:
    """ Compute the adaptive `p_response_factor`.

    Assuming that the most presented peptide
    across all simulated cells is presented e.g. 700 times on a cell, then `p_response_factor = 1/ 700`.
    This ensures that no information loss occurs.
    This method is called only if no `p_response_factor` is provided
    in the config file.
    """
    max_rep = df_cells.groupby(['cell_ids', 'presented_peptides']).size().max()
    p_response_factor = 1 / max_rep

    msg = f"adaptive `p_response_factor` is 1 / {max_rep} = {p_response_factor}"
    logger.info(msg)

    return p_response_factor
