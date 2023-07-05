""" This module contains various helpers used in the tests.

Mostly, these helpers have to do with constructing objects that will be used in
the tests.
"""
import logging
logger = logging.getLogger(__name__)

import itertools
import numpy as np
import pandas as pd

from neoag_dt import (
    AlleleScoreCache,
    Cell,
    CellFactory,
    MHC,
    Peptide,
    PeptideScoreCache,
    PeptideVaccine,
    Protein,
    ProteinCleaver,
    SomaticVariant,
    pMHC,
    TCell
)

import neoag_dt.dt_utils as dt_utils
import neoag_dt.data as dt_data

from typing import List, Sequence, Tuple, Mapping


###
# Loading from files
###
def get_hlas_df() -> pd.DataFrame():
    df_hlas = dt_data.load_hlas()
    return df_hlas


def get_peptide_sequences_df() -> pd.DataFrame:
    df_peptide_sequence = dt_data.load_peptide_sequences()
    return df_peptide_sequence


def get_proteins_df() -> pd.DataFrame:
    df_proteins = dt_data.load_proteins()
    return df_proteins


def get_variants_df() -> pd.DataFrame:
    df_variants = dt_data.load_variants()
    return df_variants


def get_simulated_cells_df() -> Tuple[pd.DataFrame, int]:
    df_cells = dt_data.load_simulated_cells()
    df_cells = df_cells[df_cells['simulation_name'] == '200_cells_10_reps']
    df_cells = df_cells[df_cells['repetition'] == 0]
    num_cells = 200
    return df_cells, num_cells


def get_distance_from_self_map() -> Mapping[str, float]:
    df_dfs = dt_data.load_distance_from_self_scores()
    dfs_map = dict(zip(
        df_dfs['Mut_peptide'],
        df_dfs['distance_from_self']
    ))
    return dfs_map


###
# Manually constructing domain objects
###
def get_protein() -> Protein:
    protein_name = 'ENSG00000197530'
    expression_mean = 7.0159899999999995
    expression_var = 1.0

    protein = Protein(
        name=protein_name,
        expression_mean=expression_mean,
        expression_var=expression_var
    )

    return protein


def get_somatic_variant(protein:Protein) -> SomaticVariant:
    # the variant is fairly common
    dna_ref_count = 5
    dna_alt_count = 10

    # there is some differential transcription in this example
    rna_ref_count = 5
    rna_alt_count = 25

    # and one of the variants from our dataset
    name = '1_1629672_C/F'

    variant = SomaticVariant(
        dna_ref_count=dna_ref_count,
        dna_alt_count=dna_alt_count,
        rna_ref_count=rna_ref_count,
        rna_alt_count=rna_alt_count,
        protein=protein,
        name=name
    )

    return variant


def get_hla_allele() -> MHC:
    name = "A0201"
    expression_mean = 57.9
    expression_var = 4.5

    protein = Protein(
        name=name, expression_mean=expression_mean, expression_var=expression_var
    )

    hla = MHC(name=name, protein=protein)
    return hla


def get_peptides() -> List[Peptide]:
    peptide_sequences = [
        'RQAEITPTK',
        'LYQIQALRW',
        'RALVVPCPR',
        'SILDNQLVR'
    ]

    peptides = [
        Peptide(peptide_sequence)
            for peptide_sequence in peptide_sequences
    ]

    return peptides


def get_hla_proteins() -> List[Protein]:
    allele_names = [
        'A0201',
        'B3501',
        'C0701',
        'A2402',
        'A0101'
    ]

    expression_mean = [
        35.2, # A0201
        0.0, # B3501, not in this sample
        5.6, # C0701
        29.9, # A2401
        0.0 # A0101, not in this sample
    ]

    expression_var = [
        2.0, # A0201
        0.1, # B3501, not in this sample
        0.5, # C0701
        3.0, # A2401
        0.1 # A0101, not in this sample
    ]

    it = zip(allele_names, expression_mean, expression_var)

    hla_proteins = [
        Protein(
            name=name, expression_mean=expression_mean, expression_var=expression_var
        ) for (name, expression_mean, expression_var) in it
    ]

    return hla_proteins


def get_hla_alleles() -> List[MHC]:
    hla_proteins = get_hla_proteins()

    hla_alleles = [
        MHC(name=p.name, protein=p) for p in hla_proteins
    ]

    return hla_alleles


def get_bound_peptide_hlas(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC]) -> List[pMHC]:

    bound_peptide_hlas = [
        phla for phla in zip(peptides, hla_alleles)
    ]

    return bound_peptide_hlas


def get_activated_t_cells(peptides_to_recognize:Sequence[Peptide]) -> List[TCell]:
    # generate a set of T cells which recognize the given peptides
    t_cells = [
        TCell(known_targets=[p], name=p.sequence)
            for p in peptides_to_recognize
    ]
    return t_cells


def get_tumor_cells(peptides_to_present:Sequence[Peptide], hla:MHC) -> List[Cell]:
    # cells should present peptide:HLA complexes with the given peptides and HLA

    # we will just create one cell for each combination of two peptides
    it = itertools.combinations(peptides_to_present, 2)

    def __get_cell(peptides):
        presented_peptides = [pMHC(p, hla) for p in peptides]
        name = ",".join([p.sequence for p in peptides])
        c = Cell(presented_peptide_hlas=presented_peptides, name=name)
        return c

    cells = [
        __get_cell(peptides) for peptides in it
    ]

    return cells


def get_optim_config(criterion: str = "MinSum") -> Mapping:
    optim_config = {
        "criterion": criterion,
        "budget": 5,
        "simulation_config": "./etc/cells-config.yaml",
        "cells_populations_file": "./analysis/cell/cell-populations.csv",
        "vaccine_element_scores": None,
        "weight_column": None,
        "peptides_file": "./neoag_dt/data/peptide-sequences.csv",
        "peptide_sequence_column": "Mut_peptide",
        "out": "./analysis/vaccines/selected-vaccine-elements.budget-5.minsum.csv",
        "p_response_factor": 0.25,
        "distance_from_self_scores": "./neoag_dt/data/distance-from-self.csv",
        "distance_from_self_column": "distance_from_self"
    }

    return optim_config


###
# Randomly constructing scorers
###
def create_random_presentation_scores(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC]) -> AlleleScoreCache:

    it = itertools.product(hla_alleles, peptides)
    rng = np.random.default_rng(seed=8675309)

    score_cache = {
        (hla.name, peptide.sequence): rng.uniform(0, 1)
            for (hla, peptide) in it
    }

    score_cache = AlleleScoreCache(score_cache, name="PresentationScoreCache")

    return score_cache


def create_random_binding_scores(
        peptides:Sequence[Peptide],
        hla_alleles:Sequence[MHC]) -> AlleleScoreCache:

    it = itertools.product(hla_alleles, peptides)
    rng = np.random.default_rng(seed=8675308)

    score_cache = {
        (hla.name, peptide.sequence): rng.uniform(3, 7)
            for (hla, peptide) in it
    }

    score_cache = AlleleScoreCache(score_cache, name="BindingScoreCache")

    return score_cache


def create_random_dfs_scores(
    peptides:Sequence[Peptide]) -> PeptideScoreCache:

    rng = np.random.default_rng(seed=8675307)

    score_cache = {
        peptide.sequence: rng.uniform(0,1)
            for peptide in peptides
    }

    score_cache = PeptideScoreCache(score_cache, name="DFSScoreCache")

    return score_cache
    

def create_random_response_scores(
    peptides:Sequence[Peptide]) -> PeptideScoreCache:

    rng = np.random.default_rng(seed=8675306)

    score_cache = {
        peptide.sequence: rng.uniform(0,1)
            for peptide in peptides
    }

    score_cache = PeptideScoreCache(score_cache, name="ResponseScoreCache")

    return score_cache


###
# Testing for vaccines
###
def get_pan_pool_1_vaccine() -> PeptideVaccine:
    peptide_sequences = [
        'AMRPNFTIK',
        'DYVYNPFMI',
        'GHFAWWTAF',
        'LPPAYTNSF',
        'RMYIFFASF',
        'VFVSNGTHW',
        'WSMATYYLF',
        'YWFFSNYLK',
    ]

    vaccine_peptides = [
        Peptide(peptide_sequence)
            for peptide_sequence in peptide_sequences
    ]

    peptide_vaccine = PeptideVaccine(
        peptides=vaccine_peptides,
        name="pan_pool_1"
    )

    return peptide_vaccine


def get_pan_pool_2_vaccine() -> PeptideVaccine:
    peptide_sequences = [
        'ALVYFLQSI',
        'ATSRTLSYY',
        'IPTNFTISV',
        'IRQEEVQEL',
        'LYLQYIRKL',
        'VFVSNGTHW',
        'YLFDESGEF',
        'YTERSEKSY',
    ]

    vaccine_peptides = [
        Peptide(peptide_sequence)
            for peptide_sequence in peptide_sequences
    ]

    peptide_vaccine = PeptideVaccine(
        peptides=vaccine_peptides,
        name="pan_pool_2"
    )

    return peptide_vaccine


def get_pan_pool_3_vaccine() -> PeptideVaccine:
    peptide_sequences = [
        'ALVYFLQSI',
        'FYGGWHNML',
        'MFVKHKHAF',
        'NFVRIIMRL',
        'RYLALYNKY',
        'SLREVRTIK',
        'VFVSNGTHW',
        'VRFPNITNL',
    ]

    vaccine_peptides = [
        Peptide(peptide_sequence)
            for peptide_sequence in peptide_sequences
    ]

    peptide_vaccine = PeptideVaccine(
        peptides=vaccine_peptides,
        name="pan_pool_3"
    )

    return peptide_vaccine


###
# Testing for BipartiteVaccineDesignModel
###
def get_ilp_model_input(criterion: str = "MinSum"):
    optim_config = get_optim_config(criterion=criterion)

    df_peptides = get_peptide_sequences_df()
    peptides = df_peptides[optim_config['peptide_sequence_column']].unique()

    df_cells, num_cells = get_simulated_cells_df()

    return optim_config, peptides, df_peptides, df_cells, num_cells, optim_config['p_response_factor']
