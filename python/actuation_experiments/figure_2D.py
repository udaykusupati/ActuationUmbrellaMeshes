import matplotlib.pyplot as plt
import numpy as np

import helpers_plots as help_
import helpers_tools as help_tools
from tools import get_center_position, height_to_percent

# ======================================================================
# ================================================================= 2D =
# ======================================================================

#------------------ GRID ------------------ 2D

def plot_undeployed_2D(input_data, curr_um,
                       show_height=False, active_cells=[], target_percents=[],
                       file_name='', show_plot=True):
    center_position = get_center_position(curr_um)

    # plot:
    fig, ax = get_ax()
    help_.ax_plot_edges(ax, input_data)
    help_.ax_annotate_index(ax, center_position)
    if show_height: help_.ax_annotate_height(ax, curr_um.umbrellaHeights, center_position)
    help_.ax_annotate_active(ax, active_cells, target_percents, center_position)

    ax.axis('equal')
    ax.axis('off')
    fig.tight_layout()
    if file_name != '':
        fig.savefig(file_name)
    plt.sca(ax)
    if show_plot: plt.show()
    else:         plt.close()

def plot2D_stress(curr_um, connectivity, init_center_pos,
                  active_cells=[], target_percents=[], stress_type='VonMises',
                  zero_as_extrem = False, show_percent=False):
    
    stress_matrix, min_, max_ = help_tools.get_smatrix_min_max(curr_um, stress_type, zero_as_extrem)
    fig, ax = get_ax()
    help_.ax_plot_stresses(ax, connectivity, stress_matrix, min_, max_, active_cells, target_percents, init_center_pos, show_percent)

    ax.axis('equal')
    ax.axis('off')
    fig.tight_layout()
    plt.sca(ax)
    plt.show()
    
def projection2D(connectivity, curr_um,
                 active_cells=[], target_percents=[], file_name='', show_plot = False):
    center_position = get_center_position(curr_um)
    fig, ax = get_ax()
    help_.ax_proj2D(ax, connectivity, active_cells, target_percents, center_position)
    
    ax.axis('equal')
    ax.axis('off')
    fig.tight_layout()
    if file_name != '':
        fig.savefig(file_name)
    plt.sca(ax)
    if show_plot: plt.show()
    else:         plt.close()

def fig_arm_stresses(connectivity,
                     active_cells,
                     percents,
                     init_center_pos,
                     show_percent,
                     s_matrix,
                     min_,
                     max_,
                     show_plot,
                     title,
                     path_names,
                     file_name=''):
    fig, ax = get_ax()
    help_.ax_plot_stresses(ax, connectivity, s_matrix, min_, max_, active_cells, percents, init_center_pos, show_percent)
    
    ax.set_title(title)
    ax.axis('equal')
    ax.axis('off')
    fig.tight_layout()
    for path in path_names:
        fig.savefig(path.format(file_name))
    plt.sca(ax)
    if show_plot: plt.show()
    else:         plt.close()
    
def fig_height2D(connectivity, active_cells, plate_thickness, init_heights, heights, positions, show_plot, title, path_names, file_name=''):
    
    indexes = list(range(len(heights))) # is this ok for regular grid ? I don't think so...
    percents = np.array(height_to_percent(init_heights, plate_thickness, indexes, heights))

    fig, ax = get_ax()
    help_.ax_dot_active_cell(ax, indexes, percents,  positions, markersize=1, edgecolor='white')
    if active_cells != []:
        help_.ax_dot_active_cell(ax, active_cells, percents[active_cells], positions, markersize=1, edgecolor='black')
    help_.ax_arms(ax, connectivity, positions)
    
    ax.set_title(title)
    ax.axis('equal')
    ax.axis('off')
    fig.tight_layout()
    for path in path_names:
        fig.savefig(path.format(file_name))
    plt.sca(ax)
    if show_plot: plt.show()
    else:         plt.close()
    


#------------------ CURVE ------------------ 1D
def figs_stress_curve(stresses, xylim, paths, show_plot=False):
    axes     = []
    path_all = []
    for path in range(len(stresses)):
        max_s = []
        ax_idx = 0
        for phase, stress_phase in enumerate(stresses[path]):
            for step, stress_step in enumerate(stress_phase):
                if phase>0 and step==0: continue
                if path==0:
                    axes.append(get_ax())
                    path_all.append([path.format(phase+1, step) for path in paths])
                
                max_s.append(max(stress_step))
                
                ax = axes[ax_idx][1]
                ax.plot(max_s, linewidth=0.85)
                plt.sca(ax)
                plt.close()
                ax_idx += 1
        
    for (fig, ax), path in zip(axes, path_all):
        ax.set_xlim(0, xylim[0])
        ax.set_ylim(0, xylim[1])
        ax.set_title('max stress per step')
        fig.tight_layout()
        for p in path:
            fig.savefig(p.format('stress_curve'))
        plt.sca(ax)
        if show_plot: plt.show()
        else:         plt.close()
    return

def figs_stress_scatter(stresses, xylim, paths, ordered=True, show_plot=False):
    ms = plt.rcParams['lines.markersize'] ** 0.5
    axes     = []
    path_all = []
    for path, stress_path in enumerate(stresses):
        ax_idx = 0
        for phase, stress_phase in enumerate(stress_path):
            for step, stress_step in enumerate(stress_phase):
                if phase>0 and step==0: continue
                if path==0:
                    axes.append(get_ax())
                    path_all.append([path.format(phase+1, step) for path in paths])
                    
                if ordered: stress_step = np.flip(np.sort(stress_step))
                
                ax = axes[ax_idx][1]
                ax.scatter(range(xylim[0]), stress_step, s=ms)
                plt.sca(ax)
                plt.close()
                ax_idx += 1
        
    for (fig, ax), path in zip(axes, path_all):
        ax.set_xlim(0, xylim[0])
        ax.set_ylim(0, xylim[1])
        if ordered: ax.set_title('max stress per arm (ordered)')
        else:       ax.set_title('max stress per arm')
        fig.tight_layout()
        for p in path:
            if ordered: fig.savefig(p.format('ordered_stress_scatter'))
            else:       fig.savefig(p.format('stress_scatter'))
        plt.sca(ax)
        if show_plot: plt.show()
        else:         plt.close()
    return

def figs_heights_curve(indexes_all, heights, xylim, active_cells_all, percents_per_steps_all, paths, show_active=True, ordered=True, show_plot=False):
    axes     = []
    path_all = []
    for path, heights_path in enumerate(heights):
        ax_idx = 0
        for phase, heights_phase in enumerate(heights_path):
            for step, heights_step in enumerate(heights_phase):
                if phase>0 and step==0: continue
                if path==0:
                    axes.append(get_ax())
                    path_all.append([path.format(phase+1, step) for path in paths])
                    
                if ordered: idx = np.flip(np.argsort(heights_step))
                else:       idx = indexes_all[path][phase]

                ax = axes[ax_idx][1]
                ax.plot(heights_step[idx], linewidth=0.85)
                if not ordered and show_active:
                    help_.ax_dot_active_cell(ax,
                                             active_cells_all[path][phase],
                                             percents_per_steps_all[path][phase][step],
                                             np.array([idx, heights_step, heights_step]).T,
                                             markersize=1)
                plt.sca(ax)
                plt.close()
                ax_idx += 1

    for (fig, ax), path in zip(axes, path_all):
        ax.set_xlim(0, xylim[0])
        ax.set_ylim(0, xylim[1])
        if ordered: ax.set_title('ordered heights per unit')
        else:       ax.set_title('heights per unit')
        fig.tight_layout()
        for p in path:
            if ordered: fig.savefig(p.format('ordered_heights_curve'))
            else:       fig.savefig(p.format('heights_curve'))
        plt.sca(ax)
        if show_plot: plt.show()
        else:         plt.close()
    return

def figs_energy_curve(energies, xylim, paths, show_plot=False): 

    axes     = []
    path_all = []
    for path, e_path in enumerate(energies):
        e = []
        ax_idx = 0
        for phase, e_phase in enumerate(e_path):
            for step, e_step in enumerate(e_phase):
                if phase>0 and step==0: continue
                if path==0:
                    axes.append(get_ax())
                    path_all.append([path.format(phase+1, step) for path in paths])
                
                e.append(e_step)
                
                ax = axes[ax_idx][1]
                ax.plot(e, linewidth=0.85)
                
                plt.sca(ax)
                plt.close()
                ax_idx += 1
                
    for (fig, ax), path in zip(axes, path_all):
        ax.set_xlim(0, xylim[0])
        ax.set_ylim(0, xylim[1])
        ax.set_title('elastic energy per step')
        fig.tight_layout()
        for p in path:
            fig.savefig(p.format('energies'))
        plt.sca(ax)
        if show_plot: plt.show()
        else:         plt.close()

def fig_empty():
    '''
    to ensure the matplotlib settings are set
    '''
    fig, ax = get_ax()
    ax.plot()
    fig.tight_layout()
    plt.sca(ax)
    plt.close()

# ---------- HELPERS ----------
def get_ax(fig_size=1, dpi=8):
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=72*dpi)
    plt.rcParams.update({'font.size': 2})
    return fig, ax