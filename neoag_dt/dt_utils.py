""" This module contains various utilities for performing the digital twin
simulations.

Many of these functions will likely be moved to more permanent homes (such as 
the `pyllars` package) in the future.
"""
import logging
logger = logging.getLogger(__name__)

import numpy as np

import typing
from typing import Optional, Tuple, Union

numpy_size = typing.Optional[typing.Union[int, typing.Tuple[int]]]
sampling_return_type = typing.Union[np.ndarray, int]

def sample_binomial(
        p:float,
        n_samples:int=1) -> sampling_return_type:
    """ Sample from a binomial distribution

     Parameters
     ----------
     p : float
        Probability of success

     n_samples : int
         The number of samples to draw from this distribution

     Returns
     -------
     samples : typing.Union[int, numpy.ndarray]
         The samples from the binomial distribution
    """
    r = np.random.binomial(n_samples, p)
    return r

def sample_beta_binomial(
        alpha:float,
        beta:float,
        n_samples:int=1,
        size:numpy_size=None) -> sampling_return_type:
    """ Sample from a beta-binomial distribution

    Parameters
    ----------
    alpha, beta : float
        The alpha and beta shape parameters of the beta distribution. They must
        be positive.

    n_samples : int
        The number of samples to draw from this distribution

    size : typing.Optional[typing.Union[int, typing.Tuple[int]]]
        The output shape of the beta distribution. Please see the
        :obj:`numpy.random.beta` documentation for details.

    Returns
    -------
    samples : typing.Union[int, numpy.ndarray]
        The samples from the beta-binomial distribution

    Note
    ----
    In order to see the outcome of multiple Bernoulli samples, the `size`
    parameter should be set to the number of desired samples, while `n_samples`
    should be set to `1`.
    """
    p = np.random.beta(alpha, beta, size=size)
    r = np.random.binomial(n_samples, p)
    return r

def sample_gamma_poisson(
        mean:float,
        var:float,
        size:numpy_size=None) -> sampling_return_type:
    """ Sample from a gamma-Poisson distribution

    Parameters
    ----------
    mean, var : float
        The respective parameters of the Gamma distribution. They must be
        non-zero, and the variance should be positive.

        **N.B.** The gamma is usually parameterized by a "shape" and "scale" 
        (usually called :math:`k` and :math:`\theta`, respectively) or a 
        "shape" and "rate" (usually called :math:`\alpha` and :math:`beta`,
        respectively). This implementation is based on simple arithmetic
        operations to derive the mean and variance from the shape and scale.

        Specifically, we have the following relationships:

        .. math:

            k = \frac{\mu^2}{\sigma^2}
            \theta = \frac{\sigma^2}{\mu}

    size : typing.Optional[typing.Union[int, typing.Tuple[int]]]
        The output shape of the beta distribution. Please see the
        :obj:`numpy.random.beta` documentation for details.

    Returns
    -------
    samples : typing.Union[int, numpy.ndarray]
        The samples from the gamma-Poisson distribution

    Note
    ----
    In order to see the outcome of multiple Poisson samples, the `size`
    parameter should be set to the number of desired samples.
    """
    shape = (mean * mean) / var
    scale = var / mean
    lam = np.random.gamma(shape, scale, size=size)
    r = np.random.poisson(lam)
    return r
