""" Create the bar chart for vaccine elements selected in the optimization

**Usage**

.. code-block::

    create-bar-chart /path/to/my/bar-charts-config.yaml --logging-level INFO

This script takes as input the selected vaccine elements from a set of
optimization settings and visualizes the frequency of how many times each
vaccine element is selected.

The script uses a yaml configuration file. The following keys are required.

* `min_hotspot_uses` : int
    The minimum number of times a vaccine element must be selected before it is
    included in the plot.

* `hue_order` : list
    The order for visualizing the populations in the plot.

`simulation_settings` : a top level dictionary. Each key is an independent
    optimization setting.  The keys of this dictionary are arbitrary and can
    be descriptive of the optimization setting.

Each optimization setting requires the following keys.

`selected_peptides` : path to optimization output
    The path to an output file from optimize-tripartite-graph-hla-population-peptides

`bar_chart` : string
    The path to the plot created by this script.
"""
import logging
logger = logging.getLogger(__name__)
import pyllars.logging_utils as logging_utils

import argparse
from typing import Mapping, Optional

# see: https://stackoverflow.com/questions/41814254
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# the order of imports seems important. Some orderings cause errors related to
# static TLS: https://github.com/scikit-learn/scikit-learn/issues/14485
import pyllars.mpl_utils as mpl_utils
import pyllars.shell_utils as shell_utils
import pyllars.utils

import pandas as pd
import seaborn as sns; sns.set(style='white', color_codes=True)

_DEFAULT_MIN_HOTSPOT_USES = 10


def create_bar_chart(
        df_selected_peptides: pd.DataFrame,
        simulation: str,
        config: Mapping,
        ax: Optional[plt.Axes]=None) -> mpl_utils.FigAx:
    """ Creates a bar chart"""

    cols = ['simulation_name', 'peptide']
    df_peptide_use = df_selected_peptides.groupby(cols).size().reset_index()
    df_peptide_use = df_peptide_use.rename(columns={0: 'count'})
    df_peptide_use = df_peptide_use.astype({'count': 'int32'})

    df_peptide_use['pretty_simulation_name'] = df_peptide_use['simulation_name']
    df_peptide_use['count'].sum()

    fig, ax = mpl_utils._get_fig_ax(ax)

    min_hotspot_uses = config.get('min_hotspot_uses', _DEFAULT_MIN_HOTSPOT_USES)

    hotspot_uses = df_peptide_use.groupby('peptide')['count'].sum()
    m_hotspot_count = hotspot_uses >= min_hotspot_uses
    common_hotspots = hotspot_uses[m_hotspot_count].index.values
    m_hotspots = df_peptide_use['peptide'].isin(common_hotspots)

    peptide_order = sorted(df_peptide_use.loc[m_hotspots, 'peptide'].unique())

    ax = sns.barplot(
        x='peptide',
        y='count',
        #x='count',
        #y='peptide',
        order=peptide_order,
        hue='pretty_simulation_name',
        data=df_peptide_use[m_hotspots],
        ax=ax,
        hue_order=config['hue_order'],
        palette=sns.color_palette("viridis", len(df_peptide_use['simulation_name'].unique()))
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    mpl_utils.set_ticklabel_rotation(ax, 60)

    legend = ax.legend(loc='upper right', title="Cells simulations")
    ax.grid(axis='y')
    ax.set_ylim((0, df_peptide_use['count'].to_numpy().max()))

    ax.set_xlabel("", fontsize=0)
    plt.xticks(rotation=90)
    ax.set_ylabel("Selection count", fontsize=20)

    sns.despine()

    return fig, ax


def _create_bar_chart(simulation, plotting_config):
    """Creates a bar chart figure for each simulation configuration and each
    population configuration"""
    msg = "creating bar charts: {}".format(simulation)
    logger.info(msg)

    simulation_config = plotting_config['simulation_settings'][simulation]
    df_selected_peptides = pd.read_csv(simulation_config['selected_peptides'])

    fig, ax = plt.subplots(figsize=(25, 5))
    fig, ax = create_bar_chart(df_selected_peptides, simulation, plotting_config, ax=ax)

    fig.suptitle(f"Vaccine optimization: {simulation}", fontsize=20)

    shell_utils.ensure_path_to_file_exists(simulation_config['bar_chart'])

    fig.savefig(simulation_config['bar_chart'], bbox_inches='tight', dpi=300)


def parse_arguments() -> argparse.Namespace:
    """Parsing function"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__
    )

    parser.add_argument('plotting_config', help="The path to the plotting "
        "configuration file")

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    return args


def main():
    """Main function to create bar charts"""
    args = parse_arguments()

    plotting_config = pyllars.utils.load_config(args.plotting_config)

    simulations = plotting_config['simulation_settings'].keys()
    for simulation in simulations:
        _create_bar_chart(simulation, plotting_config)


if __name__ == '__main__':
    main()
