import numpy as np
from matplotlib import colormaps
import matplotlib.pyplot as plt

# ======================================================================
# ============================================================== PLOTS =
# ======================================================================
def ax_arms(ax, connectivity, positions):
    arms_pos = _get_arms_pos(connectivity, positions)
    for arm in arms_pos:
        ax.plot(arm[:,0], arm[:,1], c='lightgrey', linewidth=0.5)
        
def ax_annotate_index(ax, positions):
    for i, [x,y,z] in enumerate(positions):
        ax.annotate(f'{i}', (x,y), ha='center', color='black')

def ax_annotate_height(ax, heights, positions):
    for [x,y,z], h in zip(positions, heights):
        c = _color_map(h, min(heights), max(heights))
        ax.annotate(f'[{h:.2f}]', (x,y), textcoords='offset points', xytext=(0,-2), ha='center', c=c)

def ax_annotate_active(ax, active_cells, target_percents, positions):
    for i, p in zip(active_cells, target_percents):
        r = p/100
        ax.annotate(f'{i}', positions[i][:2], ha='center', color=(r,1-r,0), weight='bold')

def ax_plot_edges(ax, input_data):
    for e in _get_edges(input_data):
        ax.plot(e[:, 0], e[:, 1], color="lightblue")

def ax_dot_active_cell(ax, active_cells, target_percents, positions, markersize=1, edgecolor='black'):
    '[RK] target_percents as color is not the wright deployment percent since phased deployement'
    'but still, this is somehow a good indication as the first phase is still ok, and each last steps too'
    s = plt.rcParams['lines.markersize']**markersize # default s value is `rcParams['lines.markersize'] ** 2`
    lw = s/20 # width of marker's edge
    for i, p in zip(active_cells, target_percents):
        r = p/100
        [x,y,_] = positions[i]
        # default zorder for plot is 2 (higher value means more on top)
        ax.scatter(x,y, color=(r,1-r,0), s=s, zorder=2.5, edgecolors=edgecolor, linewidths=lw)

def ax_plot_stresses(ax, connectivity, stress_matrix, min_, max_, active_cells, percents, positions, show_percent):
    _ax_arms_as_stress(ax, connectivity, stress_matrix, min_, max_, positions)
    ax_dot_active_cell(ax, active_cells, percents, positions, markersize=0.8)
    _ax_show_percent(ax, show_percent, active_cells, percents, positions)

def ax_proj2D(ax, connectivity, active_cells, percents, positions):
    arms_pos = _get_arms_pos(connectivity, positions)
    for arm in arms_pos:
        ax.plot(arm[:,0], arm[:,1], c='black', linewidth=0.85)
        ax_dot_active_cell(ax, active_cells, percents, positions)

# ----------------------------------------------------------------------
# ------------------------------------------------------------ helpers -
# ----------------------------------------------------------------------        
def _ax_arms_as_stress(ax, connectivity, s_matrix, min_, max_, positions):
    arms_pos = _get_arms_pos(connectivity, positions)
    colors = _get_arms_color(connectivity, s_matrix, min_, max_)
    for arm,c in zip(arms_pos, colors):
        ax.plot(arm[:,0], arm[:,1], c=c)# linewidth=8)

def _ax_show_percent(ax, show_percent, active_cells, target_percents, positions):
    for i, p in zip(active_cells*show_percent, target_percents):
        ax.annotate(f'{p: >5.0f}', positions[i][:2], ha='left', va='center', color='black')

def _get_edges(input_data):
    vertices = input_data['base_mesh_v']
    edge     = np.roll(np.insert(input_data['base_mesh_f'], 0,
                                 input_data['base_mesh_f'][:,0], axis=1),
                       -1, axis=1)
    return vertices[edge]

def _get_arms_pos(connectivity, positions):
    center_xy = np.array(list(zip(positions[:,0], positions[:,1])))
    arms_pos = center_xy[connectivity]
    return arms_pos

def _get_arms_color(connectivity, s_matrix, min_, max_):
    if   isinstance(max_, float): return [_color_map(s_matrix[i,j], min_, max_, cmap_name='') for i,j in connectivity]
    elif isinstance(max_, np.ndarray) : return [_color_map(s_matrix[i,j], 0, max_[arm]) for arm,(i,j) in enumerate(connectivity)]
    else: raise ValueError(f'max value should either be a `float` or a `np.ndarray`, is is a {type(max_)}')

def _color_map(value, min_, max_, cmap_name='', nombre=256):
    if max_ == min_:
        if max_ == 0: return 0,1,0.25      # no stress at all             -> 'green'
        else:         return 0.8, 0.8, 0.8 # max stress/height everywhere -> 'lightgrey'

    p = (value-min_)/(max_-min_)
    if cmap_name!='':
        return colormaps[cmap_name].resampled(nombre).colors[int(p*(nombre-1))]
        # do not work for 'jet' nor 'rgb, cool' or others colormap (no attribute `colors`)
    
    # default: linear red to green
    r = p #[RK] try some exponential values
    g = 1-r
    b = 0.25*g
    return r,g,b