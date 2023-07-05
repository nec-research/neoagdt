""" This module contains a vaccine element.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List, Sequence

_DEFAULT_NAME = "VaccineElement"


class VaccineElement(object):
    """ A vaccine element.

    Attributes
    ----------
    name : str
        The object name.

    weight : float
        The weight/cost of a vaccine element, considered in the ILP optimization.
    """
    def __init__(self, name: str, weight: float = 1) -> "VaccineElement":
        self.weight = weight
        self.name = name

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, s: str):
        return self.name == s

    def __ne__(self, s: str):
        return not(self == s)

    def __str__(self):
        return self.name

    @classmethod
    def construct_list(
            cls,
            element_names: Sequence[str]) -> List["VaccineElement"]:
        """ Class method to construct the object.
        """
        ve_list = [
            VaccineElement(name=n) for n in element_names
        ]

        return ve_list
