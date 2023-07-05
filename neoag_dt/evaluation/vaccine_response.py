""" Contains the function which evaluates a set of vaccines on a given population
of simulated cells.
"""
from typing import Mapping, Sequence, Tuple

import pandas as pd

from neoag_dt.evaluation.vaccine import Vaccine
from neoag_dt.evaluation.population import Population
from neoag_dt.evaluation.vaccine_evaluator import VaccineEvaluator


def evaluate_vaccines(
        simulation_name: str,
        repetition: int,
        vaccines: Sequence[Vaccine],
        population_size_map: Mapping[str, int],
        scores_dict: Mapping[Tuple[str, int], pd.DataFrame],
        df_simulations: pd.DataFrame):
    """ Evaluates a set of vaccines on a given population
    of simulated cells.
    """
    df_cells = get_simulation(df_simulations, simulation_name, repetition)
    df_scores = scores_dict[(simulation_name, repetition)]
    population_size = population_size_map[simulation_name]

    population_name = "{}, {}".format(simulation_name, repetition)
    population = Population.construct(
        df_cells, population_size, name=population_name
    )

    ret = [
        {
            "vaccine_name": v.name,
            "evaluator": VaccineEvaluator(population, v, df_scores).evaluate(progress_bar=True),
            "population_name": population_name,
            "simulation_name": simulation_name,
            "repetition": repetition
        } for v in vaccines
    ]

    return ret


def get_simulation(
        df_simulations: pd.DataFrame,
        simulation_name: str,
        repetition: int) -> pd.DataFrame:
    """ Extracts one given population simulation from the full dataframe.
    """
    # select the cells of this simulation
    simulation_groups = df_simulations.groupby(['simulation_name', 'repetition'])
    simulation_group = (simulation_name, repetition)
    df_cells = simulation_groups.get_group(simulation_group)

    return df_cells
