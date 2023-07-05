""" Calculate the likelihood of response for each cell in each population for
a set of vaccines.
Plot the statistics of the response in the investigated populations of cells.

**Usage**

.. code-block::

    evaluate-vaccine-response /path/to/my/response-likelihood-config.yaml --logging-level INFO

This script takes as input a set of vaccines, a set of populations, and
necessary scores for calculating the likelihood of response for each cell.

The script uses a yaml configuration file. The following keys are required.

* `simulations` : path to the cancer digital twins
    The path to the output from the `simulate-neoag-cancer-cells script`

* `scores` : path to scores for each vaccine element
    The path to the vaccine element scores with likelihoods of response

* `final_response_out` : path to the output file
    The file contains the likelihood of response for each cell in each
    population (that is, combination of `simulation_name` and `repetition` from
    the `simulations` file) for all of the vaccines described below.
    The evaluation considers one single final vaccine derived from each group of
    repeated optimizations.

* `sim_specific_response_out` : path to the output file
    The file contains the likelihood of response for each cell in each
    population (that is, combination of `simulation_name` and `repetition` from
    the `simulations` file) for all of the vaccines described below.
    The evaluation considers the simulation-specific vaccines from the
    repeated optimizations. This means, each population of cells is evaluated with
    the vaccine optimized for them specifically.

* `box_plot_out` : path to output figure
    The figure depicts the response distributions for the considered cell populations as
    box plots.

* `violin_plot_out` : path to output figure
    The figure depicts the response distributions for the considered cell populations as
    violin plots.

* `line_plot_out` : path to output figure
    The figure depicts the percentage of cells which respond to a vaccine ("coverage")
    as a function of a threshold.

* `sim_specific_line_plot_out`: path to output figures
    The figure depicts the percentage of cells which respond to a vaccine ("coverage")
    as a function of a threshold. It compares the final vaccine composition
    to all simulation-specific vaccines.

* `vaccines` : a dictionary specifying the vaccines to evaluate.
    Each key is a name for the respective vaccine. It is arbitrary and could be a description
    of the optimization setting. Each value is a path to a vaccine, which
    must match the output format of `optimize-vaccine-ilp`.

* `population_sizes` : a dictionary specifying the size of each simulated population.
    The keys shall be the same as the names of the simulations names.
"""
import logging

logger = logging.getLogger(__name__)
import pyllars.logging_utils as logging_utils

import argparse
import pandas as pd

from typing import List, Mapping, Tuple

import pyllars.utils
import pyllars.collection_utils as collection_utils
import pyllars.pandas_utils as pd_utils
import pyllars.shell_utils as shell_utils

from neoag_dt.evaluation.vaccine_element import VaccineElement
from neoag_dt.evaluation.vaccine import Vaccine
from neoag_dt.evaluation.vaccine_response import evaluate_vaccines
from neoag_dt.evaluation.plot import create_figures


def create_final_vaccine(vaccine_file: str, vaccine_name: str) -> Vaccine:
    """Create a vaccine from a file which only contains its vaccine elements.
    This method is meant to be used for the final vaccine compositions,
    which are derived from combining the results of multiple optimizations
    over multiple repeated simulated cells populations.
    """
    df_vaccine = pd.read_csv(vaccine_file)
    vaccine_elements = VaccineElement.construct_list(df_vaccine['peptide'])
    vaccine = Vaccine(vaccine_elements, name=vaccine_name)
    return vaccine


def create_sim_specific_vaccines(config) -> Mapping[Tuple[str, str], List[Vaccine]]:
    """For each simulation and each vaccine optimization setting, create a vaccine.
    It returns a dictionary of this type;
    <(simulation, repetition), [vaccines]>
    """
    sim_specific_vaccines = dict()

    for vaccine_name, vaccine_file in config['sim_specific_vaccines'].items():
        vaccines_df = pd.read_csv(vaccine_file)
        for simulation in vaccines_df['simulation_name'].unique():
            simulation_df = vaccines_df[vaccines_df['simulation_name'] == simulation]
            for repetition in simulation_df['repetition'].unique():
                sim_specific_vaccines[(simulation, repetition)] = []

    for vaccine_name, vaccine_file in config['sim_specific_vaccines'].items():
        vaccines_df = pd.read_csv(vaccine_file)
        for simulation in vaccines_df['simulation_name'].unique():
            simulation_df = vaccines_df[vaccines_df['simulation_name'] == simulation]
            for repetition in simulation_df['repetition'].unique():
                vaccine_df = simulation_df[simulation_df['repetition'] == repetition]
                vaccine_elements = VaccineElement.construct_list(vaccine_df['peptide'])
                name = f"{vaccine_name}.{simulation}.rep-{repetition}"
                sim_specific_vaccines[(simulation, repetition)].append(
                    (Vaccine(vaccine_elements, name=name))
                )

    return sim_specific_vaccines


def get_evaluation_results(row: pd.Series) -> pd.DataFrame:
    evaluator = row['evaluator']
    p_response = evaluator.get_p_response_values()
    log_p_response = evaluator.get_log_p_response_values()

    vaccine = evaluator.vaccine
    vaccine_name = vaccine.name

    pop = evaluator.population
    df_cells = pop.to_data_frame()

    df_cells['vaccine'] = vaccine_name
    df_cells['p_response'] = p_response
    df_cells['log_p_response'] = log_p_response

    return df_cells


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__
    )

    parser.add_argument('coverage_config', help="The configuration file for "
                                                "calculating the coverage.")

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    return args


def main():
    """Compute cell response likelihoods and evaluate vaccines efficacy.
    """
    args = parse_arguments()

    msg = "loading the coverage config file: '{}'".format(args.coverage_config)
    logger.info(msg)
    config = pyllars.utils.load_config(args.coverage_config)

    msg = "loading simulated populations: '{}'".format(config['simulations'])
    logger.info(msg)
    df_simulations = pd.read_csv(config['simulations'])
    tmp = df_simulations.groupby(['simulation_name']).nunique().reset_index('simulation_name')
    simulation_reps_map = dict(zip(tmp['simulation_name'], tmp['repetition']))

    msg = "loading population sizes"
    logger.info(msg)
    population_size_map = config['population_sizes']

    scores_dict = {}
    for sim_name in simulation_reps_map.keys():
        for rep in range(simulation_reps_map[sim_name]):
            msg = "loading vaccine element scores for sim '{}', rep '{}'".format(sim_name, rep)
            logger.info(msg)
            scores_dict[(sim_name, rep)] = pd.read_csv(config['scores'].format(sim_name, rep))

    msg = "creating the final vaccines"
    logger.info(msg)
    final_vaccines = [
        create_final_vaccine(vaccine_file, vaccine_name)
        for vaccine_name, vaccine_file in config['final_vaccines'].items()
    ]

    msg = "extracting unique populations"
    logger.info(msg)

    cols = ['simulation_name', 'repetition']
    simulation_repetitions = df_simulations[cols].drop_duplicates().values

    msg = "evaluating coverage of the final vaccines"
    logger.info(msg)

    final_eval = [
        evaluate_vaccines(
            simulation,
            repetition,
            final_vaccines,
            population_size_map,
            scores_dict,
            df_simulations
        ) for simulation, repetition in simulation_repetitions
    ]

    msg = "combining and formatting all coverage results"
    logger.info(msg)

    final_eval = collection_utils.flatten_lists(final_eval)
    df_final_eval = pd.DataFrame(final_eval)

    final_eval = pd_utils.apply(
        df_final_eval, get_evaluation_results, progress_bar=True
    )
    df_final_eval = pd.concat(final_eval)
    df_final_eval = df_final_eval.reset_index(drop=True)

    msg = "writing the results to disk: '{}'".format(config['final_response_out'])
    logger.info(msg)
    shell_utils.ensure_path_to_file_exists(config['final_response_out'])
    df_final_eval.to_csv(config['final_response_out'], index=False)

    msg = "creating the simulation-specific vaccines"
    logger.info(msg)
    sim_specific_vaccines = create_sim_specific_vaccines(config)

    msg = "evaluating coverage of the simulation-specific vaccines"
    logger.info(msg)

    sim_specific_eval = [
        evaluate_vaccines(
            simulation,
            repetition,
            sim_specific_vaccines[(simulation, repetition)],
            population_size_map,
            scores_dict,
            df_simulations
        ) for simulation, repetition in simulation_repetitions
    ]

    msg = "combining and formatting all coverage results"
    logger.info(msg)

    sim_specific_eval = collection_utils.flatten_lists(sim_specific_eval)
    df_sim_specific_eval = pd.DataFrame(sim_specific_eval)

    sim_specific_eval = pd_utils.apply(
        df_sim_specific_eval, get_evaluation_results, progress_bar=True
    )
    df_sim_specific_eval = pd.concat(sim_specific_eval)
    df_sim_specific_eval = df_sim_specific_eval.reset_index(drop=True)

    msg = "writing the results to disk: '{}'".format(config['sim_specific_response_out'])
    logger.info(msg)
    shell_utils.ensure_path_to_file_exists(config['sim_specific_response_out'])
    df_sim_specific_eval.to_csv(config['sim_specific_response_out'], index=False)

    msg = "creating figures"
    logger.info(msg)
    create_figures(df_final_eval, df_sim_specific_eval, config)


if __name__ == '__main__':
    main()
