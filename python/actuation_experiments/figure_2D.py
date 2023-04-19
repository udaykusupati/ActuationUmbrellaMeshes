import matplotlib.pyplot as plt
import numpy as np

import helpers_plots as help_
import helpers_tools as help_tools
import helpers_grid  as help_grid

FIG_SIZE_DENOM = 10

# ======================================================================
# ================================================================= 2D =
# ======================================================================

#------------------ GRID ------------------

def plot2D(input_data, curr_um,
           show_height=False, active_cells=[], target_percents=[], file_name=''):
    center_position = help_grid.get_center_position(curr_um)

    # plot:
    ax = get_ax(fig_size=curr_um.numUmbrellas()//FIG_SIZE_DENOM)
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
    ax = get_ax(fig_size=curr_um.numUmbrellas()//FIG_SIZE_DENOM)
    help_.ax_plot_stresses(ax, connectivity, stress_matrix, min_, max_, active_cells, target_percents, init_center_pos, show_percent)

    ax.axis('equal')
    plt.axis('off')
    plt.show()
    
def projection2D(connectivity, curr_um,
                 active_cells=[], target_percents=[], file_name=''):
    center_position = help_grid.get_center_position(curr_um)
    ax = get_ax(fig_size=curr_um.numUmbrellas()//FIG_SIZE_DENOM)
    help_.ax_proj2D(ax, connectivity, active_cells, target_percents, center_position)
    
    ax.axis('equal')
    plt.axis('off')
    if file_name != '':
        plt.savefig(file_name)
    plt.show()

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
    ax = get_ax(fig_size=len(init_center_pos)//FIG_SIZE_DENOM)
    help_.ax_plot_stresses(ax, connectivity, s_matrix, min_, max_, active_cells, percents, init_center_pos, show_percent)
    
    ax.set_title(title)
    ax.axis('equal')
    plt.axis('off')
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    else : plt.close()


#------------------ CURVE ------------------

def figs_heights(indexes, heights, active_cells, percents_per_steps, paths, show_active = True, show_plot=False):
    max_x, max_y, max_step, fig_size = _maxX_maxY_maxStep_figSize(heights)
    for step, heights_ in enumerate(np.array(heights).transpose((1,0,2))):
        paths_s = [path.format(step/max_step*100) for path in paths]
        
        ax = get_ax(fig_size=fig_size)
        for i,(idx, h, a) in enumerate(zip(indexes, heights_, active_cells)):
            ax.scatter(idx,h, marker='_', linewidths=plt.rcParams['lines.markersize']/FIG_SIZE_DENOM)
            if show_active: help_.ax_dot_active_cell(ax, a, percents_per_steps[i][step], np.array([idx, h, h]).T, markersize=1)
        
        ax.set_ylim(0, max_y)
        ax.set_title('heights per unit')
        for path in paths_s:
            plt.savefig(path.format('heights'))
        if show_plot: plt.show()
        else : plt.close()

def figs_stress_scatter(stresses, paths, ordered=True, show_plot=False):
    max_x, max_y, max_step, fig_size = _maxX_maxY_maxStep_figSize(stresses)
    for step, s_ in enumerate(np.array(stresses).transpose((1,0,2))):
        paths_s = [path.format(step/max_step*100) for path in paths]
        
        ax = get_ax(fig_size=fig_size)
        ms = plt.rcParams['lines.markersize'] ** 0.5
        for s in s_:
            if ordered: s = np.flip(np.sort(s))
            ax.scatter(range(len(s)), s, s=ms)
        
        ax.set_ylim(0, max_y)
        if ordered: ax.set_title('max stress per arm (ordered)')
        else : ax.set_title('max stress per arm')
        for path in paths_s:
            if ordered: plt.savefig(path.format('ordered_stress_scatter'))
            else: plt.savefig(path.format('stress_scatter'))
        if show_plot: plt.show()
        else : plt.close()

def figs_stress_curve(stresses, paths, show_plot=False):
    max_x, max_y, max_step, fig_size = _maxX_maxY_maxStep_figSize(stresses)
    max_s = []
    for step, s_ in enumerate(np.array(stresses).transpose((1,0,2))):
        paths_s = [path.format(step/max_step*100) for path in paths]
        max_s.append(np.array(s_).max(axis=1))
        
        ax = get_ax(fig_size=fig_size)
        ax.plot(max_s)
        
        ax.set_ylim(0, max_y)
        ax.set_xlim(0, max_x)
        ax.set_title('max stress per step')
        for path in paths_s:
            plt.savefig(path.format('stress_curve'))
        if show_plot: plt.show()
        else : plt.close()

def figs_energy_curve(energies, paths, show_plot=False):
    max_x, max_y, max_step, fig_size = _maxX_maxY_maxStep_figSize(energies)
    e = []
    for step, e_ in enumerate(np.array(energies).transpose()):
        paths_s = [path.format(step/max_step*100) for path in paths]
        e.append(e_)
        
        ax = get_ax(fig_size=fig_size)
        ax.plot(e)
        ax.set_ylim(0, max_y)
        ax.set_xlim(0, max_x)
        ax.set_title('elastic energy per step')
        for path in paths_s:
            plt.savefig(path.format('energies'))
        if show_plot: plt.show()
        else : plt.close()


# ---------- HELPERS ----------
def get_ax(fig_size=1, dpi=8):
    plt.rcParams.update({'font.size': 2*fig_size})
    _, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=72*dpi)
    return ax

def _maxX_maxY_maxStep_figSize(data):
    d = np.array(data)
    return d.shape[1], 1.05*d.max(), d.shape[1]-1, max(1, d.shape[0]//FIG_SIZE_DENOM)