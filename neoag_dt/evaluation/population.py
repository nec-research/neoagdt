""" This module contains a "Population" object.
"""
import logging
logger = logging.getLogger(__name__)
from typing import Sequence

import pandas as pd

from neoag_dt.evaluation.cell import Cell

_DEFAULT_NAME = "Population"


class Population(object):
    """ A population, which is primarily a set of cells.

    Attributes
    ----------
    cells : Sequence[neoag_dt.evaluation.cell.Cell]
        The cells which constitute the cancer digital twin.
    """
    def __init__(
            self,
            cells: Sequence[Cell],
            name: str = _DEFAULT_NAME) -> "Population":
        
        self.cells = cells
        self.name = name

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger
        """
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def to_data_frame(self) -> pd.DataFrame:
        """ Create a data frame from the population in which each row
        corresponds to one cell
        """
        df_cells = pd.DataFrame([
            c.to_dict() for c in self.cells
        ])

        df_cells['population'] = self.name

        return df_cells

    @classmethod
    def construct(
            cls,
            df_cells: pd.DataFrame,
            population_size: int,
            name: str) -> "Population":
        """ Class method to construct a population of cells.
        """

        msg = "building population: {}".format(name)
        logger.info(msg)
        cells = [
            Cell.construct(df_cells[df_cells['cell_ids'] == c], c)
            for c in range(population_size)
        ]

        population = Population(cells, name=name)
        return population
