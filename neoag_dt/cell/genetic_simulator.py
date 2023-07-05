""" This module contains a class to simulate the genetics for a variant.
"""
import logging
logger = logging.getLogger(__name__)

from typing import Mapping, NamedTuple

import tqdm

from neoag_dt.cell.somatic_variant import SomaticVariant
import neoag_dt.dt_utils as dt_utils


class GeneticSimulationResult(NamedTuple):
    """ The result of the genetic simulator.
    """
    variant_in_dna: bool
    protein_count: int


_DEFAULT_NAME = "GeneticSimulator"


class GeneticSimulator(object):
    """ A class to simulate the genetic processes for a variant.

    Attributes
    ----------
    {sequencing,expression}_pseudocount : float
        A pseudocount to add to all sequencing depth or expression values.These
        can be thought of as a kind of hyperprior on the parameters of the beta
        (or gamma) distribution for sampling.

        Practically, these are used to avoid passing `0`s to the sampling
        functions.

        It is straightforward to express the Jeffreys prior for the beta
        distribution (a "pseudocount" of `0.5`). It is not so straightforward to
        express the Jeffreys prior for the gamma distribution, especially if the
        observed expression is close to `0`. So we just use `1` as the default.

        Please see this discussion for more details on the Jeffreys prior for
        a gamma distribution: https://stats.stackexchange.com/questions/197173


    name : str
        A name for this object
    """
    def __init__(self,
            sequencing_pseudocount: float = 0.5,
            expression_pseudocount: float = 1,
            name: str = _DEFAULT_NAME) -> "GeneticSimulator":

        self.sequencing_pseudocount = sequencing_pseudocount
        self.expression_pseudocount = expression_pseudocount
        self.name = name

    def log(self, msg:str, level:int=logging.INFO):    
        """ Log `msg` using `level` using the module-level logger """    
        msg = "[{}] {}".format(self.name, msg)
        logger.log(level, msg)
    
    def simulate(self, variant: SomaticVariant) -> GeneticSimulationResult:
        """ Perform a genetic simulation for `variant`

        In particular, the simulation entails the following:
        * Calculate allele frequency of the variant
        * Sampling if the variant is in the DNA
        * If it is in the DNA, then sampling the protein molecule count from
            gene expression data and multiply by RNA variant allele frequency,
            since not all transcripts carry the mutation

        Parameters
        ----------
        variant : neoag_dt.SomaticVariant
            The variant. The `rna` attribute on the variant must be set.

        Returns
        -------
        simulation_result : neoag_dt.genetic_simulator.GeneticSimulationResult
            The results of the simulation. In case some steps of the simulation
            do not occur, the values are set to meaningful defaults. For
            example, if the variant is not in the RNA, then the protein count
            will always be set to 0.
        """

        protein_count = 0

        variant_in_dna = dt_utils.sample_binomial(variant.calculate_vaf_dna())
        variant_in_dna = (variant_in_dna == 1)

        if variant_in_dna:
            rna_mean = variant.protein.expression_mean + self.expression_pseudocount
            rna_var = variant.protein.expression_var
            protein_count = dt_utils.sample_gamma_poisson(rna_mean, rna_var)
            protein_count = int(protein_count * variant.calculate_vaf_rna())

        ret = GeneticSimulationResult(variant_in_dna, protein_count)
        return ret

    def simulate_all(
            self,
            variants: Mapping[str, SomaticVariant],
            progress_bar: bool = True) -> Mapping[SomaticVariant, GeneticSimulationResult]:
        """ Perform simulations for all `variants`

        Please see :obj:`~neoag_dt.genetic_simulator.GeneticSimulator.simulate`
        for details on what the simulation entails.

        Parameters
        ----------
        variants : typing.Mapping[str, neoag_dt.SomaticVariant]
            The variants

        progress_bar : bool
            Whether to show a progress bar for the simulations

        Returns
        -------
        all_simulation_results : typing.Mapping[neoag_dt.SomaticVariant, neoag_dt.genetic_simulator.GeneticSimulationResult]
            The results of the simulations
        """
        it = variants.items()
        if progress_bar:
            it = tqdm.tqdm(it)

        all_simulation_results = {
            variant : self.simulate(variant)
                for variant_name, variant in it
        }

        return all_simulation_results
