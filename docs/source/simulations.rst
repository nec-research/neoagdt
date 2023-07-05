Neoantigen Digital Twin Simulations
***********************************

This page describes the scripts to run and evaluate the Neoantigen Digital Twin
simulations. The main scripts in this package are as follows.

* `Simulating the cancer digital twin`_ : Create a simulated population of cancer cells based on the patient's sequencing data
* `Designing optimal vaccines`_ : Select vaccine elements such that the likelihood of no response is minimized for each population of cancer cells
* `Visualizing the selected vaccine elements`_ : Create a bar chart for vaccine elements selected in the optimization
* `Calculating the likelihood of response for all cells and vaccines`_ : Calculate the likelihood of response for each cancer cell in each population for a set of vaccines

Simulating the cancer digital twin
====================================

.. automodule:: neoag_dt.cli.simulate_neoag_cancer_cells
    :noindex:

Designing optimal vaccines
===========================

.. automodule:: neoag_dt.cli.optimize_vaccine_ilp
    :noindex:

Visualizing the selected vaccine elements
=========================================

.. automodule:: neoag_dt.cli.create_selected_vaccine_element_bar_chart
    :noindex:

Calculating the likelihood of response for all cells and vaccines
====================================================================

.. automodule:: neoag_dt.cli.evaluate_vaccine_response
    :noindex:
