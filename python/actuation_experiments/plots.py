import sys
sys.path.append('..')
import linkage_vis
from visualization_helper import get_color_field

import numpy as np
from matplotlib import colormaps
import matplotlib.pyplot as plt

import helpers_plots as help_
import helpers_tools as help_tools
import helpers_grid  as help_grid

FIG_SIZE = 3

# ======================================================================
# ================================================================= 3D =
# ======================================================================
def plot3D(curr_um, input_data,
           rod_colors=None, uidBased=False):
    if rod_colors is None:
        rod_colors = get_color_field(curr_um, input_data, uidBased) 
    view = linkage_vis.LinkageViewer(curr_um, width=800, height=600)
    view.update(scalarField = rod_colors)
    return view

def plot3D_stress(curr_um, stress_type,
                  verbose=True):
    view = linkage_vis.LinkageViewer(curr_um, width=800, height=600)
    stresses = help_tools.set_plate_stress_null(curr_um, help_tools.get_stresses(curr_um, stress_type))
    view.update(scalarField = stresses)
    if verbose:
        print(f'{stress_type} Stresses Extrem values:\n\
    max : {np.array(stresses).max():.2e}\n\
    min : {np.array(stresses).min():.2e}')
    return view

# ======================================================================
# ================================================================= 2D =
# ======================================================================
def plot2D(input_data, curr_um,
           show_height=False, active_cells=[], target_percents=[], file_name=''):
    center_position = help_grid.get_center_position(curr_um)

    # plot:
    ax = help_.get_ax(FIG_SIZE)
    help_.ax_plot_edges(ax, input_data)
    help_.ax_annotate_index(ax, center_position)
    if show_height: help_.ax_annotate_height(ax, curr_um.umbrellaHeights, center_position)
    help_.ax_annotate_active(ax, active_cells, target_percents, center_position)

    ax.axis('equal')
    plt.axis('off')
    # show plot
    if file_name != '':
        plt.savefig(file_name)
    plt.show()

def plot2D_stress(curr_um, connectivity, init_center_pos,
                  active_cells=[], target_percents=[], stress_type='maxBending',
                  zero_as_extrem = False, show_percent=False):
    
    stress_matrix, min_, max_ = help_tools.get_smatrix_min_max(curr_um, stress_type, zero_as_extrem)
    ax = help_.get_ax(FIG_SIZE)
    help_.ax_plot_stresses(ax, connectivity, stress_matrix, min_, max_, active_cells, target_percents, init_center_pos, show_percent)

    ax.axis('equal')
    plt.axis('off')
    plt.show()
    
def projection2D(connectivity, curr_um,
                 active_cells=[], target_percents=[], file_name=''):
    center_position = help_grid.get_center_position(curr_um)
    ax = help_.get_ax(FIG_SIZE)
    help_.ax_proj2D(ax, connectivity, active_cells, target_percents, center_position)
    
    ax.axis('equal')
    plt.axis('off')
    if file_name != '':
        plt.savefig(file_name)
    plt.show()
    