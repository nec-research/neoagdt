""" This module contains helpers for loading data files.
"""
import logging
logger = logging.getLogger(__name__)

import os
import pathlib
import yaml

import pandas as pd

from typing import AbstractSet, Mapping, Optional


def _get_base_data_dir() -> pathlib.Path:
    import neoag_dt.data # that's me!
    base_data_dir = os.path.dirname(neoag_dt.data.__file__)
    base_data_dir = pathlib.Path(base_data_dir)
    return base_data_dir

###
# Paths to files in the `neoag_dt/data` directory
###
def get_hlas_path() -> str:
    path = _get_base_data_dir() / "hlas.csv"
    return str(path)


def get_peptide_sequences_path() -> str:
    path = _get_base_data_dir() / "peptide-sequences.csv"
    return str(path)


def get_proteins_path() -> str:
    path = _get_base_data_dir() / "genes.csv"
    return str(path)


def get_variants_path() -> str:
    path = _get_base_data_dir() / "variants.csv"
    return str(path)


def get_cells_path() -> str:
    path = _get_base_data_dir() / "cell-populations.csv"
    return str(path)


def get_dfs_path() -> str:
    path = _get_base_data_dir() / "distance-from-self.csv"
    return str(path)


###
# Load data files
###
def load_hlas() -> pd.DataFrame:
    f = get_hlas_path()
    df = pd.read_csv(f)
    return df


def load_peptide_sequences(valid_lengths: Optional[AbstractSet] = {9}) -> pd.DataFrame:
    f = get_peptide_sequences_path()
    df = pd.read_csv(f)

    if valid_lengths is not None:
        msg = "Loading peptide sequences. Only keeping lengths: {}".format(
            valid_lengths
        )
        logger.info(msg)

        lengths = df['Mut_peptide'].str.len()
        m = lengths.isin(valid_lengths)
        df = df[m].reset_index(drop=True)

    return df


def load_proteins() -> pd.DataFrame:
    f = get_proteins_path()
    df = pd.read_csv(f)
    return df


def load_variants() -> pd.DataFrame:
    f = get_variants_path()
    df = pd.read_csv(f)
    return df


def load_simulated_cells() -> pd.DataFrame:
    f = get_cells_path()
    df = pd.read_csv(f)
    return df


def load_distance_from_self_scores() -> pd.DataFrame:
    f = get_dfs_path()
    df = pd.read_csv(f)
    return df
