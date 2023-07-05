""" Simulate populations of cancer cells.

**Usage**

.. code-block::

    simulate-neoag-cancer-cells  /path/to/my/cells-config.yaml --logging-level INFO

This script simulates a population of cancer cells presenting neoantigens, which
are used as input for the optimization model that selects the vaccine elements.

The script uses a yaml configuration file. The following keys are required.

* `gene file`: Path to the gene expression file. Its columns are specified by
    the keys with prefix `gene`.

* `variant_file`: Path to the variant file. Its columns are specified by the
    keys with prefix `var`.

* `peptide_file`: Path to the peptide sequence file. Its columns are specified
    by the keys with prefix `peptide`.

* `hla_file`: Path to the HLA allele file. Its columns are specified by the keys
    with prefix `hla`.

* `binding_scores_file`: Path to the file containing the predicted binding scores
    for each peptide to each HLA allele. Its columns are specified by the keys
    with prefix `binding scores`.

* `presentation_scores_file`: Path to the file containing the predicted
    presentation score for each peptide and each HLA allele.

* `cells_out`: Outputfile containing the simulated cells that will be used for the
    optimization of the vaccine.

* `gene_name_column`: str
    The name of the column that contains the Ensembl gene identifier (ENSGXY).
    Needs to be mappable to `var_gene_column`.

* `gene_expression_mean_column`: str
    The name of the column that contains the mean gene expression.

* `gene_expression_var_column`: str
    The name of the column that contains the variance gene expression. Variance
    is calculated from lower and upper bound of the 95% confidence intervall as
    given by Cufflinks.

* `gene_default_mean`: int
    The default value of the mean gene expression, should be 0.

* `gene_default_var`: int
    The default value of the gene expression variance, should be 1.

* `var_name_column`: str
    The name of the column that contains the mutation identifier consisting
    of chromosome, position and amino acid change. Needs to be mappable to
    `peptide_var_name_column`.

* `var_gene_column`: str
    The name of the column that contains the Ensembl gene identifier (ENSGXY).
    Needs to be mappable to `gene_name_column`.

* `var_dna_ref`: str
    The name of the column that contains the sequencing depth of the reference
    position in the whole exome (DNA) sequencing data.

* `var_dna_alt`: str
    The name of the column that contains the sequencing depth of the mutated
    position in the whole exome (DNA) sequencing data.

* `var_rna_ref`: str
    The name of the column that contains the sequencing depth of the reference
    position in the RNA sequencing data.

* `var_rna_alt`: str
    The name of the column that contains the sequencing depth of the mutated
    position in the RNA sequencing data.

* `vaf`: str
    The name of the column which include the DNA VAF. It can be `null`;
    in that case the VAF is computed.

* `peptide_sequence_column`: str
    The name of the column that contains the amino acid sequence of the mutated
    peptide.

* `peptide_var_name_column`: str
    The name of the column that contains the mutation identifier consisting
    of chromosome, position and amino acid change. Needs to be mappable to
    `var_name_column`.

* `hla_allele_name_column`:
    The name of the column that contains the HLA allele name in simplified
    notation (e.g. A0201)

* `hla_gene_name_column`: str
    The name of the column that contains the Ensembl gene identifier (ENSGXY)
    for each HLA gene.

* `binding_scores_allele_column_name`: str
    The name of the column that contains the HLA allele name in simplified
    notation (e.g. A0201). Needs to be one of the alleles contained in `hla_file`.

* `binding_scores_peptide_column_name`: str
    The name of the column that contains the amino acid sequence of the mutated
    peptide.

* `binding_scores_score_column_name`: str
    The name of the column that contains the binding scores (float between 0 and 1).

* `presentation_scores_allele_column_name`: str
    The name of the column that contains the HLA allele name in simplified
    notation (e.g. A0201). Needs to be one of the alleles contained in `hla_file`.

* `presentation_scores_peptide_column_name`: str
    The name of the column that contains the amino acid sequence of the mutated
    peptide.

* `presentation_scores_score_column_name`: str
    The name of the column that contains the presentation scores (float between 0 and 1).

* `simulations`: a top level dictionary
    Each key is an independent simulation setting. The keys of
    this dictionary are arbitrary and can be descriptive of the simulation setting.

Each simulation setting requires the following keys.

* `simulation_name`: str
    The name of the simulation.

* `simulation_num_cells`: int
    The number of cells that will be simulated.

* `simulation_expression_pseudocount`: int
   A pseudocount to add to all HLA expression values, default is 1.
   Further description is in the CellFactory class.

* `num_repetitions`: int
   The number of populations of cells that will be simulated.

"""
import logging
import pyllars.logging_utils as logging_utils
logger = logging.getLogger(__name__)

import argparse
import pandas as pd
import numpy as np
import random
import tqdm

import pyllars.utils
from pyllars.shell_utils import ensure_path_to_file_exists

from typing import List, Mapping, Sequence, Tuple

from neoag_dt import CellFactory
from neoag_dt.cell.protein import Protein
from neoag_dt.cell.genetic_simulator import GeneticSimulator
from neoag_dt.cell.peptide_mhc_binding_simulator import PeptideMHCBindingSimulator
from neoag_dt.cell.cell import Cell
from neoag_dt.cell.allele_score_cache import AlleleScoreCache
from neoag_dt.cell.protein_cleaver import ProteinCleaver
from neoag_dt.cell.mhc import MHC
from neoag_dt.cell.peptide import Peptide
from neoag_dt.cell.somatic_variant import SomaticVariant


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('config')

    logging_utils.add_logging_options(parser)
    parser.add_argument('--seed', type=int, default=42, help='random seed')

    args = parser.parse_args()
    logging_utils.update_logging(args)
    return args


def create_cells(
        config: Mapping,
        variant_map: Mapping[str, SomaticVariant],
        hla_alleles: Sequence[MHC],
        binding_scores: AlleleScoreCache,
        cleavage_scores: AlleleScoreCache,
        presentation_scores: AlleleScoreCache) -> List[Cell]:

    ###
    # Create the simulation objects
    ###
    msg = "Creating simulation objects"
    logger.info(msg)

    # and the object to perform the cleavage
    protein_cleaver = ProteinCleaver(presentation_scores)

    # for each variant, we need to check DNA, RNA, and protein count
    genetic_simulator = GeneticSimulator()

    # and binding competition
    binding_simulator = PeptideMHCBindingSimulator(binding_scores)

    ###
    # Create the cell actual cell factory
    ###
    cell_factory = CellFactory(
        genetic_simulator=genetic_simulator,
        protein_cleaver=protein_cleaver,
        binding_simulator=binding_simulator,
        presentation_scores=presentation_scores,
        expression_pseudocount=config.get('simulation_expression_pseudocount'),
        name='cell_factory'
    )

    ###
    # Finally create the cells
    ###
    msg = "Simulating the cells"
    logger.info(msg)

    cells = [
        cell_factory.create_cell(variant_map, hla_alleles, progress_bar=False)
            for _ in tqdm.trange(config.get('simulation_num_cells'))
    ]

    return cells


def simulation_to_df(
        populations: List[List[Cell]],
        peptide_map: Mapping[str, Peptide],
        simulation_name: str) -> pd.DataFrame:
    repetitions = []
    cell_ids = []
    presented_peptides = []
    presented_hlas = []

    for repetition, population in enumerate(populations):
        for cell_id, cell in enumerate(population):
            presented_peptides += [
                pep_hla.peptide.sequence
                for pep_hla in cell.presented_peptide_hlas
            ]
            presented_hlas += [
                pep_hla.hla.name
                for pep_hla in cell.presented_peptide_hlas
            ]

            cell_ids += [cell_id]*len(cell.presented_peptide_hlas)
            repetitions += [repetition]*len(cell.presented_peptide_hlas)

    d = {
        'repetition': repetitions,
        'cell_ids': cell_ids,
        'presented_peptides': presented_peptides,
        'presented_hlas': presented_hlas,
        'simulation_name': [simulation_name]*len(repetitions)
    }
    df = pd.DataFrame(data=d)
    df['mutation'] = df['presented_peptides'].apply(lambda x: peptide_map[x].somatic_variant)
    return df


def intersect_mutations(
        df_var: pd.DataFrame, df_pep: pd.DataFrame, config: Mapping
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mut_var = set(df_var[config.get('var_name_column')].unique())
    mut_pep = set(df_pep[config.get('var_name_column')].unique())

    intersection = mut_var.intersection(mut_pep)
    xor = mut_var ^ mut_pep

    logger.info(f"Mutations {xor} excluded because not common to both variants AND peptides dataframes")

    logger.info(f"Mutations {intersection} considered because common to both variants AND peptides dataframes")

    df_var = df_var[df_var[config.get('var_name_column')].isin(intersection)]
    df_pep = df_pep[df_pep[config.get('var_name_column')].isin(intersection)]

    return df_var, df_pep


def main():
    args = parse_args()
    config = pyllars.utils.load_config(args.config)

    np.random.seed(args.seed)
    random.seed(args.seed)

    ###
    # File IO
    ###
    df_genes = pd.read_csv(config['gene_file'])
    df_variants = pd.read_csv(config['variant_file'])
    df_peptide_sequences = pd.read_csv(config['peptide_file'])
    df_hlas = pd.read_csv(config['hla_file'])
    df_binding_scores = pd.read_csv(config['binding_scores_file'])
    df_presentation_scores = pd.read_csv(config['presentation_scores_file'])

    # consider only mutations present in both peptides and variants files
    df_variants, df_peptide_sequences = intersect_mutations(
        df_variants, df_peptide_sequences, config
    )

    ###
    # Create the domain objects
    ###    
    msg = "Creating domain objects"
    logger.info(msg)

    protein_map = Protein.create_proteins(
        df_genes,
        name_column=config.get('gene_name_column'),
        expression_mean_column=config.get('gene_expression_mean_column'),
        expression_var_column=config.get('gene_expression_var_column'),
        default_mean_value=config.get('gene_default_mean'),
        default_var_value=config.get('gene_default_var')
    )

    variant_map = SomaticVariant.create_somatic_variants(
        df_variants,
        protein_map=protein_map,
        name_column=config.get('var_name_column'),
        gene_column=config.get('var_gene_column'),
        dna_ref_count_column=config.get('var_dna_ref'),
        dna_alt_count_column=config.get('var_dna_alt'),
        rna_ref_count_column=config.get('var_rna_ref'),
        rna_alt_count_column=config.get('var_rna_alt'),
        vaf_column=config.get('vaf')
    )

    peptide_map = Peptide.create_peptides(
        df_peptide_sequences,
        variant_map=variant_map,
        sequence_column=config.get('peptide_sequence_column'),
        variant_column=config.get('peptide_var_name_column'),
    )

    hla_alleles = MHC.create_hlas(
        df_hlas,
        protein_map=protein_map,
        allele_name_column=config.get('hla_allele_name_column'),
        gene_name_column=config.get('hla_gene_name_column'),
    )

    ###
    # Build the score caches
    ###
    msg = "Building the binding and presentation score caches"
    logger.info(msg)
    
    binding_scores = AlleleScoreCache.construct(
        df_binding_scores,
        allele_column=config.get('binding_scores_allele_column_name'),
        peptide_column=config.get('binding_scores_peptide_column_name'),
        score_column=config.get('binding_scores_score_column_name'),
        name="binding_scores"
    )

    cleavage_scores = AlleleScoreCache.construct(
        df_presentation_scores,
        allele_column=config.get('cleavage_scores_allele_column_name'),
        peptide_column=config.get('cleavage_scores_peptide_column_name'),
        score_column=config.get('cleavage_scores_score_column_name'),
        name="cleavage_scores"
    )
    
    presentation_scores = AlleleScoreCache.construct(
        df_presentation_scores,
        allele_column=config.get('presentation_scores_allele_column_name'),
        peptide_column=config.get('presentation_scores_peptide_column_name'),
        score_column=config.get('presentation_scores_score_column_name'),
        name="presentation_scores"
    )

    simulations_dfs = []
    for simulation in config['simulations']:
        ###
        # Run the simulation
        ###
        msg = f"Starting cancer cells simulations `{simulation['simulation_name']}`: " \
              f"{simulation['simulation_num_cells']} cells x {simulation['num_repetitions']} repetitions"
        logger.info(msg)

        populations = [
            create_cells(
                config=simulation,
                variant_map=variant_map,
                hla_alleles=hla_alleles,
                binding_scores=binding_scores,
                cleavage_scores=cleavage_scores,
                presentation_scores=presentation_scores
            )
            for _ in range(simulation['num_repetitions'])
        ]
        simulations_dfs.append(simulation_to_df(populations, peptide_map, simulation['simulation_name']))

    ###
    # Save the simulation results
    ###
    msg = f"Writing simulated cells to disk: '{config['cells_out']}'"
    logger.info(msg)
    df = pd.concat(simulations_dfs)
    ensure_path_to_file_exists(config['cells_out'])
    df.to_csv(config['cells_out'], index=False)


if __name__ == '__main__':
    main()
