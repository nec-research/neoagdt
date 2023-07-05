"""This module contains a "Cell" object.
This is a more lightweight version compared to `neoag_dt.cell.cell.Cell`.
"""
import logging
logger = logging.getLogger(__name__)
from typing import List, Dict, Union

import pandas as pd

_DEFAULT_NAME = "Cell"


class Cell(object):
    """ A cell,
    which is primarily a set of presented peptides.

    Attributes
    ----------
    presented_peptides : List[str]
        The list of presented peptides.

    num_presented_peptides : int
        The number of presented peptides.

    name : str
        Object name.
    """
    def __init__(
            self,
            presented_peptides: List[str],
            num_presented_peptides: int,
            name: int = _DEFAULT_NAME) -> "Population":
        
        self.presented_peptides = presented_peptides
        self.num_presented_peptides = num_presented_peptides
        self.name = name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, cell_name: int):
        return self.name == cell_name

    def __ne__(self, cell_name: int):
        return not (self == cell_name)

    def __int__(self):
        return self.name

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger
        """
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def to_dict(self) -> Dict:
        """ Create a dictionary representing this cell.
        """
        ret = {
            'name': self.name,
            'num_presented_peptides': self.num_presented_peptides
        }

        return ret

    @classmethod
    def construct(
            cls,
            df_presented_pMHCs: pd.DataFrame,
            name: Union[str, int]) -> "Population":
        """ Class method to construct a cell."""

        presented_peptides = df_presented_pMHCs['presented_peptides'].unique().tolist()
        num_presented_peptides = len(presented_peptides)
        cell = Cell(presented_peptides, num_presented_peptides, name=name)
        return cell
