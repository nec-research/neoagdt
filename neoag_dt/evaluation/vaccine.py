""" This module contains a class for a "Vaccine" domain object.
"""
import logging
logger = logging.getLogger(__name__)
from typing import Sequence

from neoag_dt.evaluation.vaccine_element import VaccineElement

_DEFAULT_NAME = "Vaccine"


class Vaccine(Sequence):
    """ A vaccine, which is primarily a set of vaccine elements.

    Attributes
    ----------
    vaccine_elements : Sequence[neoag_dt.evaluation.vaccine_element.VaccineElement]
        The vaccine elements constituting the vaccine.

    name : str
        The object name.
    """
    def __init__(
            self,
            vaccine_elements: Sequence[VaccineElement],
            name: str = _DEFAULT_NAME) -> "Vaccine":
        
        self.vaccine_elements = vaccine_elements
        self.name = name
        super().__init__()

    def __getitem__(self, i: int):
        return self.vaccine_elements[i]

    def __len__(self):
        return len(self.vaccine_elements)

    def log(self, msg: str, level: int = logging.INFO):
        """ Log `msg` using `level` using the module-level logger
        """
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)
