import numpy as np
from matplotlib import colormaps
import matplotlib.pyplot as plt

import helpers_grid as help_grid

# ======================================================================
# ============================================================== PLOTS =
# ======================================================================
def ax_plot_edges(ax, input_data):
    for e in _get_edges(input_data):
        ax.plot(e[:, 0], e[:, 1], color="lightblue")

def ax_annotate_index(ax, center_position):
    for i, [x,y,z] in enumerate(center_position):
        ax.annotate(f'{i}', (x,y), ha='center', color='black')

def ax_annotate_height(ax, curr_um, center_position):
    um_heights = curr_um.umbrellaHeights
    for [x,y,z], h in zip(center_position, um_heights):
        c = _color_map(h, min(um_heights), max(um_heights))
        ax.annotate(f'[{h:.2f}]', (x,y), textcoords='offset points', xytext=(0,-12), ha='center', c=c)

def ax_annotate_active(ax, active_cells, target_percents, center_position):
    for i, p in zip(active_cells, target_percents):
        r = p/100
        ax.annotate(f'{i}', center_position[i][:2], ha='center', color=(r,1-r,0), weight='bold')

def ax_plot_stresses(ax, connectivity, stress_matrix, min_, max_, active_cells, percents, position, show_percent):
    _ax_arms_as_stress(ax, connectivity, stress_matrix, min_, max_, position)
    _ax_dot_active_cell(ax, active_cells, percents, position)
    _ax_show_percent(ax, show_percent, active_cells, percents, position)

def ax_proj2D(ax, input_data, active_cells, percents, position):
    arms_pos = _get_arms_pos(input_data, position)
    for arm in arms_pos:
        ax.plot(arm[:,0], arm[:,1], c='black')#, linewidth=8)
        _ax_dot_active_cell(ax, active_cells, percents, position)

def fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent, s_matrix, min_, max_, show_plot, title, path_names, file_name=''):
    ax = get_ax()
    ax_plot_stresses(ax, connectivity, s_matrix, min_, max_, active_cells, percents, init_center_pos, show_percent)
    ax.set_title(title)
    ax.axis('equal')
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()

def fig_stress_curve(max_stresses, max_x, max_y, show_plot, title, path_names, file_name=''):
    ax = get_ax()
    ax.plot(max_stresses)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)
    ax.set_title(title)
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()

def fig_stress_scatter(s_matrix, max_y, show_plot, title, path_names, file_name=''):
    ax = get_ax()
    sort = np.flip(np.unique(s_matrix[s_matrix!=0]))
    ax.scatter(range(len(sort)), sort, s=5)
    ax.set_ylim(0, max_y)
    ax.set_title(title)
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()

def get_ax(fig_size=3, dpi=8*72):
    _, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=dpi)
    if fig_size == 3: # is_default
        plt.rcParams.update({'font.size': 5})
    return ax

# ----------------------------------------------------------------------
# ------------------------------------------------------------ helpers -
# ----------------------------------------------------------------------
def _ax_dot_active_cell(ax, active_cells, target_percents, positions):
    s = plt.rcParams['lines.markersize'] ** 1.4 # default s value is `rcParams['lines.markersize'] ** 2`
    lw = s/15 # width of marker's edge
    for i, p in zip(active_cells, target_percents):
        r = p/100
        [x,y,_] = positions[i]
        ax.scatter(x,y, color=(r,1-r,0), s=s, zorder=2.5, edgecolors='white', linewidths=lw) # default zorder for plot is 2 (higher means more on top)
        

def _ax_arms_as_stress(ax, connectivity, s_matrix, min_, max_, position):
    arms_pos = _get_arms_pos(connectivity, position)
    colors = _get_arms_color(connectivity, s_matrix, min_, max_)
    for arm,c in zip(arms_pos, colors):
        ax.plot(arm[:,0], arm[:,1], c=c)#, linewidth=8)

def _ax_show_percent(ax, show_percent, active_cells, target_percents, position):
    for i, p in zip(active_cells*show_percent, target_percents):
        ax.annotate(f'{p: >5.0f}', position[i][:2], ha='left', va='center', color='black')

def _get_edges(input_data):
    vertices = input_data['base_mesh_v']
    edge     = np.roll(np.insert(input_data['base_mesh_f'], 0,
                                 input_data['base_mesh_f'][:,0], axis=1),
                       -1, axis=1)
    return vertices[edge]

def _get_arms_pos(connectivity, position):
    center_xy = np.array(list(zip(position[:,0], position[:,1])))
    arms_pos = center_xy[connectivity]
    return arms_pos

def _get_arms_color(connectivity, s_matrix, min_, max_):
    try: # chek if max_ is float or not:
        int(max_)
        return [_color_map(s_matrix[i,j], min_, max_, cmap_name='') for i,j in connectivity]
    except TypeError: # it is not
        return [_color_map(s_matrix[i,j], 0, max_[arm]) for arm,(i,j) in enumerate(connectivity)]

def _color_map(value, min_, max_, cmap_name='', nombre=256):
    if max_ == min_:
        if max_ == 0: return 0,1,0.25 # no stress at all -> "green"
        else:         return 0,0,0    # max stress/height everywhere

    p = (value-min_)/(max_-min_)
    if cmap_name!='':
        return colormaps[cmap_name].resampled(nombre).colors[int(p*(nombre-1))]
        # do not work for 'jet' nor 'rgb, cool' or others colormap (no attribute `colors`)
    
    # default: linear red to green
    r = p #[RK] try some exponential values
    g = 1-r
    b = 0.25*g
    return r,g,b