Running a simple example
========================

This package includes the necessary configuration and data files to run a
simple example of the digital twin analysis. The example consists of four main
steps:

1. `Simulate the cancer digital twin`_
2. `Optimize vaccine designs according to the simulated cancer digital twin`_
3. `Evaluate the efficacy of the vaccines on the cancer digital twin`_
4. `Plot the commonly-selected vaccine elements`_

This document describes how to perform each analysis step, as well as the
expected input and output for each step. All of the commands below assume the package has been successfully installed
and that execution is from the base folder of the package.

**N.B.** Several of the input files used for the examples below are **randomly
generated**. That is, the files are syntactically correct, but they contain
nonsense values. **They should not be used for actual analysis.**

Simulate the cancer digital twin
------------------------------------

The first step of the analysis entails simulation of a set of cancer cells, i.e.
the cancer digital twin. The ``simulate-cancer-cells`` script performs these
simulations.

The repository contains sample input files in ``neoag_dt/data``.
The major input files for this script are summarized in the following:

* genes file (``/neoag_dt/data/genes.csv``): it contains gene expression information.
* variants file (``/neoag_dt/data/variants.csv``): it contains information about the patient's variants.
* peptides file (``/neoag_dt/data/peptide-sequences.csv``): it contains the peptides associated to each variant.
* HLA file (``/neoag_dt/data/hlas.csv``): the patient's HLA alleles.
* binding scores file (``/neoag_dt/data/peptide-hla-scores.csv``): predicted binding scores for each peptide to each HLA allele.
* presentation scores file (``/neoag_dt/data/peptide-hla-scores.csv``): file containing the predicted presentation score for each peptide and each HLA allele.

The paths to these files can be set in a dedicated configuration file, which is
passed via command line (see below). Please see the :doc:`./simulations` documentation for an exact
specification of this file. The repository includes the ``cells-config.yaml`` file, which can be
used as template.

The cancer cell simulation simulates transcription, translation, intracellular processing,
HLA binding and cell surface presentation.

.. code-block::

    simulate-cancer-cells etc/cells-config.yaml --logging-level INFO

After running this program, the code will create the ``analysis/cell/cell-populations.csv``
file. It will contain 10 simulated populations, of 1000 cancer cells each (if the
config file has not been altered), where each digital twin is a population of cancer cells.
Each cell is a set of presented peptide-MHC (pMHC) complexes.

Optimize vaccine designs according to the simulated cancer digital twin
--------------------------------------------------------------------------

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
proxy (``/neoag_dt/data/distance-from-self.csv`` or leave as ``None`` to skip).
Afterwards, the optimization selects the optimal subset of vaccine elements which maximizes
the likelihood of immune response based on the approximated scores. The optimization consists
in a max flow problem, solved via integer linear programming (ILP).

.. code-block::

    optimize-vaccine-ilp etc/optimization-config.yaml --num-procs <NUM_PROCS> --num-threads-per-proc <NUM_THREADS_PER_PROC> --logging-level INFO

``--num-procs <NUM_PROCS>``: number of parallel processes (and of CPU cores).
``--num-threads-per-proc <NUM_THREADS_PER_PROC>``: number of threads to allocate for each process.
For example: ``--num-procs 10 --num-threads-per-proc 1``. This way, you use 10 CPU cores, solve
10 optimization problems in parallel, 1 problem per CPU. It best to set ``--num-threads-per-proc 1``.

The script creates the files ``analysis/vaccines/selected-vaccine-elements.budget-30.minsum.csv``
and ``analysis/vaccines/vaccine.budget-30.minsum.csv``. The first file includes
the optimal vaccine elements for each cancer cell population. The second file
contains the final vaccine elements, pooled from the various optimization based on
how many times they were selected. These files give the optimal vaccine designs for the
respective optimization settings.

Evaluate the efficacy of the vaccines on the cancer digital twin
---------------------------------------------------------------------

The package includes a script for evaluating the expected coverage likelihood
for each cancer digital twin.

.. code-block::

    evaluate-vaccine-response etc/response-likelihood-config.yaml --logging-level INFO

This script creates the ``analysis/evaluation/sim-specific-response-likelihoods.csv``
and ``analysis/evaluation/final-response-likelihoods.csv`` files.

The first file contains the likelihood of response for each cell in each
population. The evaluation considers the simulation-specific vaccines from the
repeated optimizations. This means, each population of cells is evaluated with
the vaccine optimized for them specifically.

The second file contains the likelihood of response for each cell in each population, just like the first one.
However, the evaluation considers one single final vaccine derived from each group of
repeated optimizations, i.e. the final pooled vaccine.

Related figures are stored in ``analysis/evaluation/plots``.

Plot the commonly-selected vaccine elements
-------------------------------------------

The package includes a script which characterizes the frequency of selected
vaccine elements. Specifically, it counts and visualizes as a bar chart the
vaccine elements which are included in a vaccine for each population more than
some minimum number of times.

.. code-block::

    create-bar-chart etc/bar-charts-config.yaml --logging-level INFO

The script creates the ``analysis/vaccines/bar-chart.budget-30.minsum.png`` figure.