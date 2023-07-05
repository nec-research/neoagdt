""" Utilities to create figures and plot for the vaccine evaluation.
"""
import logging
logger = logging.getLogger(__name__)

from typing import List

import pyllars.shell_utils as shell_utils
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

sns.set(style='white', color_codes=True)
sns.set_palette("magma", 1)


def create_figures(df_final_eval: pd.DataFrame, df_sim_specific_eval: pd.DataFrame, config):
    """ Creates figures.
    """
    vaccines = df_final_eval['vaccine'].unique().tolist()

    df_plot = df_final_eval.copy(deep=True)
    df_plot = df_plot.drop(
        columns=['num_presented_peptides', 'p_response', 'log_p_response', 'vaccine']
    )
    df_plot.drop_duplicates(inplace=True)
    df_plot.drop(columns=['name'], inplace=True)
    df_plot['population'] = df_plot['population'].apply(lambda x: x.split(", ")[0])
    df_plot.sort_values(by=['population'], inplace=True)

    for v in vaccines:
        df_results_vaccine = df_final_eval.copy(deep=True)
        df_results_vaccine = df_results_vaccine[df_results_vaccine['vaccine'] == v]
        df_results_vaccine = df_results_vaccine.reset_index()
        df_results_vaccine['population'] = df_results_vaccine['population'].apply(lambda x: x.split(", ")[0])
        df_results_vaccine.sort_values(by=['population'], inplace=True)
        df_plot[v] = df_results_vaccine['p_response']

    df_melted = _melt_df(df_plot, vaccines)
    _create_bar_plots(df_melted, config)
    _create_violin_plots(df_melted, config)
    _create_box_plots(df_melted, config)
    _plot_coverage_curves(df_final_eval, config)
    _plot_sim_specific_coverage_curves(df_final_eval, df_sim_specific_eval, config)


def _melt_df(df: pd.DataFrame, vaccines: List[str]):
    """ Melts the dataframe with the response probabilities for
    the each (population, vaccine).
    """
    df_melted = pd.melt(
        df,
        id_vars=['population'],
        value_vars=vaccines,
        var_name='vaccines',
        value_name='Cell response probability'
    )
    df_melted = df_melted.rename(columns={'population': 'Cells simulations'})
    return df_melted


def _create_bar_plots(df: pd.DataFrame, config) -> None:
    """ Creates bar plots for a set of vaccines and populations.
    """
    fig, ax = plt.subplots()
    ax = sns.barplot(
        x='Cells simulations',
        y='Cell response probability',
        data=df,
        hue='vaccines',
        ci='sd',
        ax=ax,
        palette=sns.color_palette("magma", len(config['final_vaccines'].keys()))
    )
    ax.legend(
        # bbox_to_anchor=(0, 1.02, 1, 0.2),
        # loc="lower left",
        loc="best",
        # mode="expand",
        title='Vaccine optimization',
        # ncol=len(config['final_vaccines'].keys())
    )
    ax.grid(axis='y')
    plt.tight_layout()
    plt.ylim(0, 1)
    shell_utils.ensure_path_to_file_exists(config['bar_plot_out'])
    plt.savefig(config['bar_plot_out'], dpi=300)


def _create_violin_plots(df: pd.DataFrame, config) -> None:
    """ Creates violin plots for a set of vaccines and populations.
    """
    fig, ax = plt.subplots()
    ax = sns.violinplot(
        x='Cells simulations',
        y='Cell response probability',
        data=df,
        hue='vaccines',
        ax=ax,
        palette=sns.color_palette("magma", len(config['final_vaccines'].keys()))
    )
    ax.legend(
        # bbox_to_anchor=(0, 1.02, 1, 0.2),
        # loc="lower left",
        loc="best",
        # mode="expand",
        title='Vaccine optimization',
        # ncol=len(config['final_vaccines'].keys())
    )
    plt.tight_layout()
    shell_utils.ensure_path_to_file_exists(config['violin_plot_out'])
    plt.savefig(config['violin_plot_out'], dpi=300)


def _create_box_plots(df: pd.DataFrame, config) -> None:
    """ Creates box plots for a set of vaccines and populations.
    """
    fig, ax = plt.subplots()
    sns.boxplot(
        x='Cells simulations',
        y='Cell response probability',
        data=df,
        hue='vaccines',
        ax=ax,
        palette=sns.color_palette("magma", len(config['final_vaccines'].keys()))
    )
    ax.legend(
        # bbox_to_anchor=(0, 1.02, 1, 0.2),
        # loc="lower left",
        loc="best",
        # mode="expand",
        title='Vaccine optimization',
        # ncol=len(config['final_vaccines'].keys())
    )
    plt.tight_layout()
    shell_utils.ensure_path_to_file_exists(config['box_plot_out'])
    plt.savefig(config['box_plot_out'], dpi=300)


def _plot_coverage_curves(df: pd.DataFrame, config) -> None:
    """ Creates coverage curves (one for each final vaccine composition), plotting
    the percentage of population which responds to the vaccine
    for different thresholds over the probability of response.
    """
    fig, ax = plt.subplots()
    thresholds = np.linspace(0, 1, num=100)
    vaccines = df['vaccine'].unique().tolist()
    vaccine_coverage_map = dict()

    for v in vaccines:
        vaccine_coverage_map[v] = _coverage_by_threshold(
            thresholds, df, v
        )

    df_plot = pd.DataFrame(vaccine_coverage_map)
    df_plot['Threshold'] = thresholds
    df_plot = df_plot.set_index('Threshold')
    df_plot = df_plot.rename_axis('Vaccine optimization', axis=1)
    sns.lineplot(
        data=df_plot,
        ax=ax,
        palette=sns.color_palette("magma", len(config['final_vaccines'].keys())),
        linewidth=3
    )
    ax.grid(axis='y')
    ax.grid(axis='x')
    ax.set_ylabel("Coverage ratio")
    shell_utils.ensure_path_to_file_exists(config['line_plot_out'])
    plt.savefig(config['line_plot_out'], dpi=300)


def _plot_sim_specific_coverage_curves(
        df_final_eval: pd.DataFrame, df_sim_specific_eval: pd.DataFrame, config
) -> None:
    """ Creates coverage curves (one for each combination of
    vaccine, simulation, repetition). It plots
    the percentage of population which responds to the
    simulation-specific vaccine compared to the final vaccine
    for different thresholds over the probability of response.
    """
    thresholds = np.linspace(0, 1, num=100)
    final_vaccines = df_final_eval['vaccine'].unique().tolist()
    sim_specific_vaccines = df_sim_specific_eval['vaccine'].unique().tolist()
    for f_v in final_vaccines:
        for sim in config['population_sizes'].keys():
            fig, ax = plt.subplots()
            vaccine_coverage_map = {'Threshold': [], 'Coverage ratio': []}

            for s_s_v in sim_specific_vaccines:
                if f_v in s_s_v:
                    vaccine_coverage_map['Coverage ratio'] += _coverage_by_threshold(
                        thresholds, df_sim_specific_eval, s_s_v
                    )
                    vaccine_coverage_map['Threshold'] += thresholds.tolist()

            df_sim_specific = pd.DataFrame(vaccine_coverage_map)
            ax = sns.lineplot(
                data=df_sim_specific,
                #palette=['k'],
                x='Threshold',
                y='Coverage ratio',
                ax=ax
            )

            vaccine_coverage_map = dict()
            vaccine_coverage_map['Final vaccine composition'] = _coverage_by_threshold(
                thresholds, df_final_eval, f_v
            )
            df_final = pd.DataFrame(vaccine_coverage_map)
            df_final['Threshold'] = thresholds
            df_final = df_final.set_index('Threshold')
            df_final = df_final.rename_axis('Vaccine', axis=1)
            ax = sns.lineplot(
                data=df_final,
                palette=['black'],
                ax=ax
            )
            ax.grid(axis='y')
            ax.grid(axis='x')
            # where some data has already been plotted to ax
            handles, labels = ax.get_legend_handles_labels()
            patch = Line2D([0], [0], label='Repeated optimizations (95% CI)', color=sns.color_palette("magma", 1)[0])
            handles.append(patch)
            plt.legend(handles=handles, loc='upper right')

            num_cells = config['population_sizes'][sim]
            reps = _get_reps(df_sim_specific_eval, sim)
            title = f"{num_cells} cells | {reps} repetitions \n Optimization: {f_v}"
            ax.set_title(title, y=1.04)
            ax.set_ylabel("Coverage ratio")

            save_path = config['sim_specific_line_plot_out'].format(f_v, sim)
            shell_utils.ensure_path_to_file_exists(save_path)
            plt.savefig(save_path, dpi=300)


def _coverage_by_threshold(
        thresholds: np.ndarray, df: pd.DataFrame, vaccine: str
) -> List[float]:
    """ Compute the percentage of responding cells for different
    threshold over the response likelihood.
    """
    coverage = []
    df_vaccine = df[df['vaccine'] == vaccine]
    for t in thresholds:
        p = len(df_vaccine[df_vaccine['p_response'] > t]) / len(df_vaccine)
        coverage.append(p)
    return coverage


def _get_reps(df_sim_specific: pd.DataFrame, sim: str) -> int:
    """ Get number of repeated simulations.
    """
    reps = 0
    for simulation in df_sim_specific['population'].unique():
        if sim in simulation:
            reps += 1
    return reps
