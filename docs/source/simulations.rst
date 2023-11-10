Neoantigen Digital Twin Simulations
***********************************

This page describes the scripts to run and evaluate the Neoantigen Digital Twin
simulations including explanations of the configuration files . The main scripts in this package are as follows.

* `Simulating the cancer digital twin`_ : Create a simulated population of cancer cells based on the patient's sequencing data
* `Designing optimal vaccines`_ : Select vaccine elements such that the likelihood of no response is minimized for each population of cancer cells
* `Visualizing the selected vaccine elements`_ : Create a bar chart for vaccine elements selected in the optimization
* `Calculating the likelihood of response for all cells and vaccines`_ : Calculate the likelihood of response for each cancer cell in each population for a set of vaccines

Simulating the cancer digital twin
====================================

.. automodule:: neoag_dt.cli.simulate_neoag_cancer_cells
    :noindex:

Based on the configurations in the ``cells-config.yaml`` file, the cells will be simulated taking different steps of the antigen presentation pathway into account. For each cell, presence of a mutation if so its expression, cleavage, binding to MHC and presentation on the cell surface is calculated by distribution sampling, producing a set of presented peptides that compose one cell.
The config file requires the following variables to be set:

Paths to files:

* ``gene_file``: Path to the file that contains gene expression information inlucding a gene identifier that is shared with the variant file.
* ``variant_file``: Path to the file that contains variant information either as VAF or as read counts
* ``peptide_file``: Path to the file that contains peptides sequences
* ``hla_file``: Path to the file that contains the HLA information and the gene identifier for each allele
* ``binding_scores_file``: Path to the file containing binding prediction scores
* ``presentation_scores_file``: Path to the file containing presentation scores. This can be the same as the binding prediction scores file as columns containing the respective information are specified later.
* ``cells_out``: Path to the CSV file that contains the cell population as a list of presented peptides and their associated cell identifier

Columns of ``gene_file``:

* ``gene_name_column``: Column containing gene identifiers
* ``gene_expression_mean_column``: Column containing mean gene expression
* ``gene_expression_var_column``: Column containing gene expression variance
* ``gene_default_mean``: Column containing default mean gene expression
* ``gene_default_var``: Column containing default gene expression variance

Columns of ``variant_file``:

* ``var_name_column``: Column containing mutation identifiers
* ``var_gene_column``: Column containing gene identifiers
* ``var_dna_ref``: Column containing tumor DNA reference read count can be set to null if VAF is given
* ``var_dna_alt``: Column containing tumor DNA mutation read count can be set to null if VAF is given
* ``var_rna_ref``: Column containing tumor RNA reference read count can be set to null if VAF is given
* ``var_rna_alt``: Column containing tumor RNA mutation read count can be set to null if VAF is given
* ``vaf``: Column containing variant allele frequency can be set to null if read counts are given

Columns of ``peptide_file``:

* ``peptide_sequence_column``: Column containing mutated peptides
* ``peptide_var_name_column``: Column containing mutation identifiers

Columns of ``hla_file``:

* ``hla_allele_name_column``: Column containing the HLA alleles
* ``hla_gene_name_column``: Column containing the gene identifier used in ``gene_file``

Columns of ``binding_scores_file``:

* ``binding_scores_allele_column_name``: Column containing the HLA alleles
* ``binding_scores_peptide_column_name``: Column containing mutated peptides
* ``binding_scores_score_column_name``: Column containing the MHC binding prediction scores
* ``cleavage_scores_allele_column_name``: Column containing the HLA alleles
* ``cleavage_scores_peptide_column_name``: Column containing mutated peptides
* ``cleavage_scores_score_column_name``: Column containing the cleavage prediction scores

Columns of ``presentation_scores_file``:

* ``presentation_scores_allele_column_name``: Column containing the HLA alleles
* ``presentation_scores_peptide_column_name``: Column containing mutated peptides
* ``presentation_scores_score_column_name``: Column containing the pMHC presentation scores

Configuration of the cell simulations:

* ``simulations``: For each simulation, the following variables need to be defined:
* ``simulation_name``: Unique name or identifier of the simulation
* ``simulation_num_cells``: Number of cells to be simulated
* ``simulation_expression_pseudocount``: Expression pseudocount value
* ``num_repetitions``: Number of simulation repetitions

Designing optimal vaccines
===========================

.. automodule:: neoag_dt.cli.optimize_vaccine_ilp
    :noindex:

The second step of the analysis is meant to select the optimal vaccine elements for each
of the simulated populations. In addition to the population, the optimization algorithm also includes a "budget".
By default, the budget indicates the number of vaccine elements which can be included.

The optimization algorithm is, generally, deterministic. Thus, it is not
normally needed to run it multiple times for the same population, set of scores,
and budget.

The repository contains the configuration file ``optimization-config.yaml``.
Most interestingly, the optimization can be solved in the peptides space, setting ``vaccine_elements: peptides``,
or in the mutation space, setting ``vaccine_elements: mutations``.
Since multiple populations of cells were simulated at the previous step, multiple
optimization problems are solved, one for each population.

First, the optimization module estimates the likelihood of immune response for the vaccine
elements, either peptides, or mutations. The core intuition is that, the more
a peptide is presented, the higher is the likelihood that including that peptide (or the
relative mutation) in a vaccine will elicit immune response. Optionally, the optimization module
can also approximate TCR recognition, by using distance-from-self as a
proxy (``distance-from-self.csv`` or leave as ``None`` to skip).
Afterwards, the optimization selects the optimal subset of vaccine elements which maximizes
the likelihood of immune response based on the approximated scores. The optimization consists
in a max flow problem, solved via integer linear programming (ILP).

The config file ``optimization-config.yaml`` contains the following variables to be set:

* ``optimization_settings``: For each optimization, the following variables need to be defined:
* ``optimization_name``: This is the name of the optimzation that serves as its unique downstream identifier
* ``criterion``: Optimization critetion, has to be either ``MinSum`` or ``MinMax``
* ``budget``: Number of elements that should be contained in the final vaccine
* ``max_solving_time``: Maximum of number of seconds spent for solving
* ``simulation_config``: Path to the simulation config file
* ``cells_populations_file``: Path to the cell populations file
* ``approximate_vaccine_element_scores``: Path to the file where the approximate vaccine element scores are stored. Use ``{}`` as placeholder for population size and number of repetition
* ``distance_from_self_scores``: Path to file containing scores for distance to self (if available)
* ``distance_from_self_column``: Column name containing scores for distance to self (if available)
* ``weight_column``: Column name containing weights for distance to self (if available)
* ``peptides_file``: Path to the file that contains peptides sequences
* ``peptide_sequence_column``: Column containing mutated peptides
* ``out``: Path to output file containing the selected vaccine elements from each simulation
* ``vaccine_out``: Path to output file containing the actual vaccine composition
* ``p_response_factor``: Response factor, can be set to null
* ``vaccine_elements``: Can be either ``mutations`` or ``peptides``
* ``mutation_id_column``: Column containing mutation identifiers

Calculating the likelihood of response for all cells and vaccines
====================================================================

.. automodule:: neoag_dt.cli.evaluate_vaccine_response
    :noindex:

This package evaluates the vaccine response for a given vaccine composition using the simulated cell populations. The configurations are stored in ``response-likelihood-config.yaml``, which requires the following varibles to be set:

* ``simulations``: Path to the cell populations file
* ``scores``: Path to the file where the approximate vaccine element scores are stored. Use ``{}`` as placeholder for population size and number of repetition
* ``final_response_out``: Path to the output containing final response likelihoods
* ``sim_specific_response_out``: Path to the output containing simulation specific response likelihoods
* ``box_plot_out``: Path to the box plot output png
* ``violin_plot_out``: Path to the violin plot output png
* ``bar_plot_out``: Path to the bar plot output png
* ``line_plot_out``: Path to line plot output png
* ``sim_specific_line_plot_out``: Path to the simulation specific line plot output pngs. Use ``{}`` as placeholder for population size and number of repetition
* ``final_vaccines``: For each vaccine optimization, the following variable needs to be set:
* ``optimization_name``: Path to output file containing the actual vaccine composition
* ``sim_specific_vaccines``: For each vaccine optimization, the following variable needs to be set:
* `optimization_name``: Path to output file containing the selected vaccine elements from each simulation
* ``population_sizes``: For each population size of the simulation, the following variable needs to be set:
* ``simulation_name`` from cells-config.yaml: Number of cells

Visualizing the selected vaccine elements
=========================================

.. automodule:: neoag_dt.cli.create_selected_vaccine_element_bar_chart
    :noindex:

The package includes a script which characterizes the frequency of selected vaccine elements. Specifically, it counts and visualizes as a bar chart the vaccine elements which are included in a vaccine for each population more than some minimum number of times. The configuration file ``bar-charts-config.yaml`` contains the following variables to be set:

* ``min_hotspot_uses``: Number of hotspot uses
* ``hue_order``: List indicated by dashes with elements consisting of each simulation
* ``simulation_settings``: For each simulation, the following variables need to be set:
* ``simulation_name``: Name of the simulation
* ``selected_peptides``: Path to file containing the selected vaccine elements from each simulation
* ``bar_chart``: Path to bar chart output png