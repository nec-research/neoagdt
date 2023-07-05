""" Select vaccine elements such that the likelihood of no response is minimized
for each population

**Usage**

.. code-block::

    optimize-vaccine-ilp  /path/to/my/optimization-config.yaml --num-procs <NUM_PROCS> --num-threads-per-proc <NUM_THREADS_PER_PROC> --logging-level INFO

This script takes as input the simulated digital twin populations of cells, the
initial set of epitopes and the calculated likelihoods of no response.
The optimization selects the vaccine elements which optimize our coverage criteria.

The script uses a yaml configuration file. The following keys are required.

* `optimization_settings` : a top level dictionary. Each key is an independent
    optimization setting. The keys of this dictionary are arbitrary and can
    be descriptive of the optimization setting.

Each optimization setting requires the following keys.

* `budget` : float
    The budget for the optimization. By default, each vaccine element has a
    cost of 1. This can be adapted programmatically but not through the
    command line interface.

* `optimization_config` : int
    Max time available (seconds) for ILP optimization.

* `criterion` : str
    The criterion (objective function) to optimize. We support `MinSum` and `MinMax`.
    `MinSum`: selects vaccine elements in order to maximize the likelihood to kill
    all cancer cells (in our ILP implementation, we minimize the likelihood of having
    no response, that's why it's called `MinSum`).
    `MinMax`: selects vaccine elements in order to minimize the likelihood of no response
    for the cell which has the highest no response likelihood.

* `cells_populations_file` : path to the digital twin cancer cells
    The path to the output from the simulate-cancer-cells command.

* `simulation_config`: path to the simulation config file
    The path to the file used to obtain `cells_simulation`

* `peptide_file` : path to the input peptides
    The path to the csv file which contains the input peptides.

* `peptide_sequence_column` : str
    The name of the column of the `peptide_file` which contains the input peptides.

* `approximate_vaccine_element_scores` : path where the approximate score files are written
    The approximated P(R=- | v_i, c_j) are stored here.

* `weight_column` : str
    The name of the column which contains the
    weight assigned to each peptides. If it is `null`, we assume all peptides
    have a constant weight of 1.

* `out` : string
    The path to the selected vaccine elements for this optimization setting.

* `vaccine_out` : string
    The path to the final vaccine composition.

* `p_response_factor` : float or None
    If a cancer cell c_j presents a peptide v_i N times, including vaccine element v_i
    in a vaccine guarantees an immune response defined by:
    `P(R=+ | v_i, c_j) = min(1, N * p_response_factor) - epsilon`
    If this parameter is `None` (`null` in the yaml config file), then the system computes an adaptive `p_response_factor`
    from the simulation output. More specifically, assuming that the most presented peptide
    across all simulated cells is presented e.g. 700 times on a cell, then `p_response_factor = 1/ 700`.
    This ensures that no information loss occurs.

* `distance_from_self_scores` : path to distance-from-self score
    This csv file shall contain for each neoantigen a distance-from-self score.
    These scores are used by the optimization as a proxy of the TCR recognition.
    These scores are multiplicative factors which contribute to the estimation
    of P(R=+ | v_i, c_j).
    The expected scores shall be in [0, 1) range. If the path is `none`,
    these distance-from-self scores are not considered.
    `P(R=+ | v_i, c_j) = (min(1, N * p_response_factor) - epsilon) * distance_from_self(v_i)`

* `distance_from_self_column` : string
    The name of the column with distance-from-self scores in the
    `distance_from_self_scores` csv file.

* `vaccine_elements` : string
    If this is "peptides", run optimization over minimal epitopes.
    If it is "mutations", run optimization over mutations.

* `mutation_id_column` : string
    The column name of the mutations.
    It can be None is `vaccine_elements = peptides`.
"""
import logging

logger = logging.getLogger(__name__)
import pyllars.logging_utils as logging_utils
from typing import Mapping, Optional, Tuple, Union

import argparse
import numpy as np
import pandas as pd
import random
from timeit import default_timer as timer

import pyllars.dask_utils as dask_utils
import pyllars.pandas_utils as pd_utils
import pyllars.shell_utils as shell_utils
import pyllars.utils

from neoag_dt.optimization.bipartite_ilp_model import BipartiteVaccineDesignModel
from neoag_dt.optimization import optimization_utils


def get_weights(
        df_peptides: pd.DataFrame,
        peptide_sequence_column: str,
        weight_column: Optional[str]) -> Mapping[str, float]:
    peptides_weights_map = None

    if weight_column is None:
        peptides_weights_map = {
            p: 1 for p in df_peptides[peptide_sequence_column].unique()
        }

    else:
        peptides_weights_map = pd_utils.dataframe_to_dict(
            df_peptides,
            key_field=peptide_sequence_column,
            value_field=weight_column
        )

    return peptides_weights_map


def select_peptides(
        simulation_group: Tuple[str, int, int],
        optimization_config: Mapping[str, Union[str, float, int]],
        logging_args: Optional[Mapping[str, str]] = None) -> pd.DataFrame:
    if logging_args is not None:
        logging_utils.update_logging(logging_args)

    msg = "optimizing simulation: {}".format(simulation_group)
    logger.info(msg)

    msg = "loading simulated_cells"
    logger.info(msg)
    df_simulations = pd.read_csv(optimization_config['cells_populations_file'])
    df_simulations = df_simulations[df_simulations['simulation_name'] == simulation_group[0]]
    df_cells = df_simulations[df_simulations['repetition'] == simulation_group[1]]

    msg = "loading vaccine element candidates"
    logger.info(msg)
    df_peptides = pd.read_csv(optimization_config['peptides_file'])
    peptides = df_peptides[optimization_config['peptide_sequence_column']].unique()
    mutations = df_peptides[optimization_config['mutation_id_column']].unique()

    if optimization_config['vaccine_elements'] == 'peptides':
        vaccine_elements = peptides
        vaccine_element_column = optimization_config['peptide_sequence_column']
    else:
        vaccine_elements = mutations
        vaccine_element_column = optimization_config['mutation_id_column']

    msg = "approximating log-probabilities of no response logP(R=- | v_i, c_j)"
    logger.info(msg)
    df_scores = optimization_utils.approximate_response_scores(
        df_cells,
        peptides,
        simulation_group[2],
        optimization_config
    )

    if optimization_config['vaccine_elements'] == "mutations":
        df_scores = optimization_utils.aggregate_by_mutations(
            df_scores, df_peptides, optimization_config
        )

    out = optimization_config['approximate_vaccine_element_scores']
    out = out.format(simulation_group[0], simulation_group[1])
    shell_utils.ensure_path_to_file_exists(out)
    df_scores.to_csv(out, index=False)

    msg = "creating weight map for each vaccine element"
    logger.info(msg)
    weight_column = optimization_config.get('weight_column')
    weights_map = get_weights(
        df_peptides, vaccine_element_column, weight_column
    )

    ilp_model = BipartiteVaccineDesignModel(
        config=optimization_config,
        peptides=vaccine_elements,
        df_cells=df_cells,
        df_scores=df_scores,
        peptide_weights_map=weights_map,
        vaccine_elements=optimization_config['vaccine_elements'],
        seed=logging_args.seed,
    )

    ilp_model.build_model()
    start = timer()
    ilp_model.optimize(max_solving_time=optimization_config['max_solving_time'])
    elapsed = timer() - start

    msg = "checking optimization results"
    logger.info(msg)
    ilp_model.check_optimization()

    msg = "retrieving selected vaccine elements"
    logger.info(msg)
    df_selected_peptides = ilp_model.get_selected_peptides()

    # add weight
    df_selected_peptides['weight'] = df_selected_peptides['peptide'].apply(
        lambda x: weights_map[x]
    )

    # add elapsed time
    df_selected_peptides['run_time'] = elapsed

    return df_selected_peptides


def run_optimization(
        simulation_name: str,
        simulation_settings: Mapping,
        args: argparse.Namespace):
    msg = "performing simulation: {}".format(simulation_name)
    logger.info(msg)

    msg = "loading the simulation config file"
    logger.info(msg)
    simulation_config = simulation_settings['simulation_config']
    simulation_config = pyllars.utils.load_config(simulation_config)

    # create the list of keys for all simulation groups
    simulation_group_vals = optimization_utils.get_simulation_group_vals(simulation_config)

    msg = "connecting to the dask cluster"
    logger.info(msg)
    dask_client, cluster = dask_utils.connect(args)

    msg = "running jobs on the dask cluster"
    logger.info(msg)
    all_peptide_dfs = dask_utils.apply_iter(
        simulation_group_vals,
        dask_client,
        select_peptides,
        simulation_settings,
        logging_args=args,
        return_futures=False,
        progress_bar=True
    )

    dask_client.close()

    msg = "combining selected data frames"
    logger.info(msg)
    df_selected_peptides = pd.concat(all_peptide_dfs)
    df_selected_peptides = df_selected_peptides.reset_index(drop=True)

    msg = "writing selected vaccine elements to disk"
    logger.info(msg)

    out = simulation_settings['out']
    shell_utils.ensure_path_to_file_exists(out)
    df_selected_peptides.to_csv(out, index=False)

    df_vaccine = optimization_utils.extract_best_peptides(
        df_selected_peptides,
        simulation_settings['budget']
    )
    out = simulation_settings['vaccine_out']
    shell_utils.ensure_path_to_file_exists(out)
    df_vaccine.to_csv(out, index=False)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__
    )

    parser.add_argument('simulations', help="The path to the simulation "
                                            "configuration file")

    dask_utils.add_dask_options(parser)
    logging_utils.add_logging_options(parser)
    parser.add_argument('--seed', type=int, default=42, help='random seed')

    args = parser.parse_args()
    logging_utils.update_logging(args)
    return args


def main():
    args = parse_arguments()

    np.random.seed(args.seed)
    random.seed(args.seed)

    msg = "loading the optimization config file"
    logger.info(msg)
    simulation_config = pyllars.utils.load_config(args.simulations)

    simulations = simulation_config['optimization_settings']
    for optimization_name, optimization_settings in simulations.items():
        start = timer()
        run_optimization(optimization_name, optimization_settings, args)
        msg = f"Elapsed time for ({optimization_name}, {optimization_settings}): {timer()-start}"
        logger.info(msg)


if __name__ == '__main__':
    main()
