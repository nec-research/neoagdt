""" A set of tests for the digital twin idea
"""
import pytest

import logging
logger = logging.getLogger(__name__)

import pyllars.logging_utils as logging_utils
logging_utils.set_logging_values(logging_level='DEBUG')

import itertools
import numpy as np
import pandas as pd

from neoag_dt import (
    AlleleScoreCache,
    Cell,
    CellFactory,
    MHC,
    Peptide,
    Protein,
    ProteinCleaver,
    SomaticVariant
)

import neoag_dt.dt_utils as dt_utils
import neoag_dt.data as dt_data
import neoag_dt.test_utils as test_utils

from typing import List, Sequence, Tuple
pHLA = Tuple[Peptide, MHC]
    
@pytest.fixture
def protein() -> Protein:
    return test_utils.get_protein()

@pytest.fixture
def somatic_variant(protein:Protein) -> SomaticVariant:
    return test_utils.get_somatic_variant(protein)

@pytest.fixture
def hla_allele() -> MHC:
    return test_utils.get_hla_allele()

@pytest.fixture
def hla_alleles() -> Sequence[MHC]:
    return test_utils.get_hla_alleles()

@pytest.fixture
def peptides() -> Sequence[Peptide]:
    return test_utils.get_peptides()

@pytest.fixture
def presentation_scores(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC]) -> AlleleScoreCache:
    return test_utils.create_random_presentation_scores(peptides, hla_alleles)

@pytest.fixture
def binding_scores(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC]) -> AlleleScoreCache:
    return test_utils.create_random_binding_scores(peptides, hla_alleles)

@pytest.fixture
def bound_peptide_hlas(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC]) -> Sequence[pHLA]:
    return test_utils.get_bound_peptide_hlas(peptides, hla_alleles)


def test_binomial() -> None:
    p = 0.6
    n_samples = 1

    sample = dt_utils.sample_binomial(
        p=p, n_samples=n_samples
    )

    assert isinstance(sample, int)

def test_beta_binomial() -> None:
    alpha = 20
    beta = 25
    n_samples = 1

    sample = dt_utils.sample_beta_binomial(
        alpha=alpha, beta=beta, size=n_samples
    )
    
    n_samples = 20
    many_samples = dt_utils.sample_beta_binomial(
        alpha=alpha, beta=beta, size=n_samples
    )

    assert isinstance(sample, np.ndarray)
    assert isinstance(many_samples, np.ndarray)

def test_gamma_poisson_reparameterization() -> None:
    mean = 4.5
    var = 2
    size = (20)

    samples = dt_utils.sample_gamma_poisson(mean, var, size)
    assert isinstance(samples, np.ndarray)
    

def test_somatic_variants() -> None:
    dna_ref_count = 5
    dna_alt_count = 10
    rna_ref_count = 20
    rna_alt_count = 25
    protein = None

    variant = SomaticVariant(
        dna_ref_count=dna_ref_count,
        dna_alt_count=dna_alt_count,
        rna_ref_count=rna_ref_count,
        rna_alt_count=rna_alt_count,
        protein=protein
    )

    assert(variant.rna_alt_count == 25)

def test_somatic_variant_sampling(somatic_variant:SomaticVariant) -> None:

    ###
    # select variants in DNA
    ###
    # BetaBernoulli(counts)
    dna_alpha = somatic_variant.dna_ref_count
    dna_beta = somatic_variant.dna_alt_count
    variant_in_dna = dt_utils.sample_beta_binomial(dna_alpha, dna_beta)
    
    # if the variant is not in the DNA, we should stop here

    ###
    #  select variants in RNA
    ###
    # BetaBernoulli(count)
    rna_alpha = somatic_variant.rna_ref_count
    rna_beta = somatic_variant.rna_alt_count
    variant_in_rna = dt_utils.sample_beta_binomial(rna_alpha, rna_beta)
    
    # if the variant is not in the RNA, we should stop here

    ###
    # determine abundance of protein
    ###
    # gamma-Poisson(TPM)
    protein_mean = somatic_variant.protein.expression_mean
    protein_var = somatic_variant.protein.expression_var
    protein_count = dt_utils.sample_gamma_poisson(protein_mean, protein_var)

def test_hla_sampling(hla_allele:MHC) -> None:
    ###
    # determine abundance of HLA molecules
    ###
    # gamma-Poisson(TPM)
    hla_mean = hla_allele.protein.expression_mean
    hla_var = hla_allele.protein.expression_var
    hla_count = dt_utils.sample_gamma_poisson(hla_mean, hla_var)

    assert isinstance(hla_count, int)

def test_protein_cleavage(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC],
        presentation_scores:AlleleScoreCache,
        num_proteins:int=20) -> None:

    # for each peptide
    #   find its maximum likelihood of presentation for all relevant alleles
    #
    # these are the weights for sampling
    peptide_weights = np.array([
        presentation_scores.get_maximum_peptide_score(peptide, hla_alleles)
            for peptide in peptides
    ])

    # normalize the sampling distribution to sum to 1
    total_weight = np.sum(peptide_weights)
    peptide_weights = peptide_weights / total_weight

    # for each protein molecule
    #   sample the peptide based on the weighted 
    selected_peptides = np.random.choice(
        peptides,
        size=num_proteins,
        replace=True,
        p=peptide_weights
    )

    assert len(selected_peptides) == num_proteins

def test_protein_cleaver(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC],
        presentation_scores:AlleleScoreCache,
        num_proteins:int=20) -> None:

    protein_cleaver = ProteinCleaver(presentation_scores)
    selected_peptides = protein_cleaver.cleave(
        peptides, hla_alleles, num_proteins
    )

    assert len(selected_peptides) == num_proteins

def test_hla_peptide_binding(
        peptides:Sequence[Peptide],
        hla_allele:MHC,
        binding_scores:AlleleScoreCache,
        num_alleles:int=20) -> None:

    # for each peptide
    #   find its likelihood of binding to this allele
    #
    # these weights will be used for sampling
    peptide_weights = np.array([
        binding_scores.get_score(hla_allele, peptide)
            for peptide in peptides
    ])

    # normalize the sampling distribution to sum to 1
    total_weight = np.sum(peptide_weights)
    peptide_weights = peptide_weights / total_weight

    # for each allele
    #   sample the peptide to which it binds based on the weighted 
    selected_peptides = np.random.choice(
        peptides,
        size=num_alleles,
        replace=True,
        p=peptide_weights
    )

    assert len(selected_peptides) == num_alleles

def test_peptide_hla_presentation(
        bound_peptide_hlas:Sequence[pHLA],
        presentation_scores:AlleleScoreCache) -> None:

    # for each bound complex
    #   check if the complex is presented
    presentation_likelihoods = np.array([
        presentation_scores.get_score(hla, peptide)
            for peptide, hla in bound_peptide_hlas
    ])

    presented = np.random.binomial(n=1, p=presentation_likelihoods)
    num_presented = np.sum(presented)

    # mostly a dummy test
    assert(num_presented <= len(bound_peptide_hlas))


def main():

    protein = test_utils.get_protein()
    somatic_variant = test_utils.get_somatic_variant(protein)    
    hla_allele = test_utils.get_hla_allele()
    
    hla_alleles = test_utils.get_hla_alleles()
    peptides = test_utils.get_peptides()
    
    presentation_scores = test_utils.create_random_presentation_scores(peptides, hla_alleles)
    binding_scores = test_utils.create_random_binding_scores(peptides, hla_alleles)    
    bound_peptide_hlas = test_utils.get_bound_peptide_hlas(peptides, hla_alleles)

    test_binomial()
    test_beta_binomial()
    test_gamma_poisson_reparameterization()
    test_somatic_variants()

    test_somatic_variant_sampling(somatic_variant)

    test_hla_sampling(hla_allele)
    
    test_protein_cleavage(peptides, hla_alleles, presentation_scores)
    test_protein_cleaver(peptides, hla_alleles, presentation_scores)

    test_hla_peptide_binding(peptides, hla_allele, binding_scores)

    test_peptide_hla_presentation(bound_peptide_hlas, presentation_scores)

if __name__ == '__main__':
    main()