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
    if show_height: help_.ax_annotate_height(ax, curr_um, center_position)
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
    
# ======================================================================
# ======================================================== SAVE IMAGES =
# ======================================================================
def plot1D_steps(connectivity, active_cells, percents_per_steps, init_center_pos, stresses_per_steps, dir_name='./00_test',
                 stress_type='maxBending', show_percent=False, show_plot=False):
    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()
    max_x, max_y = steps, 1.05*max_stress_all
    
    title = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
    path_names = []
    if dir_name == '': img_format = []
    else:              img_format = ['jpg', 'png']
    for f in img_format:
        path_names.append(f'{dir_name}/{f}/{stress_type}_'+'{{}}_{:0>3.0f}Deployment'+f'.{f}')
    
    # random perturbations does affect undeployed state
    deployed = False
    
    for s, (s_matrix, percents) in enumerate(zip(stresses_per_steps,percents_per_steps)):
        max_stress_step = s_matrix[s_matrix != 0].max()

        title_s = title.format(s/steps*100, max_stress_step)
        path_names_s = []
        for path in path_names:
            path_names_s.append(path.format(s/steps*100))
        
        help_.fig_stress_scatter(connectivity, deployed*s_matrix, max_y, show_plot, title_s.format('stress'), path_names_s, 'scatter', ordered=False)
        
        # perturbations do not affect deployed state
        deployed = True

def plot2D_steps(connectivity, active_cells, percents_per_steps, init_center_pos, stresses_per_steps, dir_name='./00_test',
                 stress_type='maxBending', show_percent=False, show_plot=False):
    
    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()
    max_x, max_y = steps, 1.05*max_stress_all

    src = np.array(connectivity)[:,0]
    dst = np.array(connectivity)[:,1]
    max_stress_per_arm = stresses_per_steps.transpose()[src,dst].max(axis=1)
    
    title = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
    path_names = []
    if dir_name == '': img_format = []
    else:              img_format = ['jpg', 'png']
    for f in img_format:
        path_names.append(f'{dir_name}/{f}/{stress_type}_'+'{{}}_{:0>3.0f}Deployment'+f'.{f}')
    
    # random perturbations does affect undeployed state
    deployed = False
    
    for s, (s_matrix, percents) in enumerate(zip(stresses_per_steps,percents_per_steps)):
        min_stress_step = s_matrix[s_matrix != 0].min()
        max_stress_step = s_matrix[s_matrix != 0].max()
        max_stresses.append(max_stress_step)

        title_s = title.format(s/steps*100, max_stress_step)
        path_names_s = []
        for path in path_names:
            path_names_s.append(path.format(s/steps*100))
        
        # normalized with general extrems values
        help_.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_all, deployed*max_stress_all, show_plot, title_s.format('all'), path_names_s, 'structure_all')
        
        # normalized with step extrems values
        help_.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_step, deployed*max_stress_step, show_plot, title_s.format('step'), path_names_s, 'structure_perSteps')

        # normalized with own extrems values
        help_.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, 0, deployed*max_stress_per_arm, show_plot, title_s.format('own'), path_names_s, 'structure_own')
        
        # ordered stresses
        help_.fig_stress_scatter(connectivity, deployed*s_matrix, max_y, show_plot, title_s.format('ordered stress'), path_names_s, 'ordered_scatter')
        
        # stress curve
        help_.fig_stress_curve(deployed*max_stresses, max_x, max_y, show_plot, title_s.format('deployment'), path_names_s, 'sPlot')
        
        # perturbations do not affect deployed state
        deployed = True

