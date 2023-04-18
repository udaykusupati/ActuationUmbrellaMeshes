import numpy as np
from matplotlib import colormaps
import matplotlib.pyplot as plt

# ======================================================================
# ============================================================== PLOTS =
# ======================================================================
def ax_annotate_index(ax, positions):
    for i, [x,y,z] in enumerate(positions):
        ax.annotate(f'{i}', (x,y), ha='center', color='black')

def ax_annotate_height(ax, heights, positions):
    for [x,y,z], h in zip(positions, heights):
        c = _color_map(h, min(heights), max(heights))
        ax.annotate(f'[{h:.2f}]', (x,y), textcoords='offset points', xytext=(0,-12), ha='center', c=c)

def ax_annotate_active(ax, active_cells, target_percents, positions):
    for i, p in zip(active_cells, target_percents):
        r = p/100
        ax.annotate(f'{i}', positions[i][:2], ha='center', color=(r,1-r,0), weight='bold')

def ax_plot_edges(ax, input_data):
    for e in _get_edges(input_data):
        ax.plot(e[:, 0], e[:, 1], color="lightblue")

def ax_plot_stresses(ax, connectivity, stress_matrix, min_, max_, active_cells, percents, positions, show_percent):
    _ax_arms_as_stress(ax, connectivity, stress_matrix, min_, max_, positions)
    _ax_dot_active_cell(ax, active_cells, percents, positions)
    _ax_show_percent(ax, show_percent, active_cells, percents, positions)

def ax_proj2D(ax, connectivity, active_cells, percents, positions):
    arms_pos = _get_arms_pos(connectivity, positions)
    for arm in arms_pos:
        ax.plot(arm[:,0], arm[:,1], c='black')#, linewidth=8)
        _ax_dot_active_cell(ax, active_cells, percents, positions)

def fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent, s_matrix, min_, max_, show_plot, title, path_names, file_name=''):
    ax = get_ax()
    ax_plot_stresses(ax, connectivity, s_matrix, min_, max_, active_cells, percents, init_center_pos, show_percent)
    ax.set_title(title)
    ax.axis('equal')
    plt.axis('off')
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    else : plt.close()


def figs_heights(indexes, heights, active_cells, percents_per_steps, paths, show_active = True, show_plot=False):
    max_y = 1.05*np.array(heights).max()
    steps = np.array(heights).shape[1]-1
    for s, (heights_, p_) in enumerate(zip(np.array(heights).transpose((1,0,2)), np.array(percents_per_steps).transpose((1,0,2)))):
        paths_s = [path.format(s/steps*100) for path in paths]
        ax = get_ax()
        ax.set_ylim(0, max_y)
        for idx, h, a, p in zip(indexes, heights_, active_cells, p_):
            ax.scatter(idx,h, marker='_')
            if show_active: _ax_dot_active_cell(ax, a, p, np.array([idx, h, h]).T)
        # ax.set_title(title)
        for path in paths_s:
            plt.savefig(path.format('heights'))
        if show_plot: plt.show()
        else : plt.close()


def figs_stress_curve(stresses, paths, show_plot=False):
    max_y = 1.05*np.array(stresses).max()
    max_x = np.array(stresses).shape[1]
    steps = max_x-1
    max_s = []
    for step, s_ in enumerate(np.array(stresses).transpose((1,0,2))):
        paths_s = [path.format(step/steps*100) for path in paths]
        max_s.append(np.array(s_).max(axis=1))
        ax = get_ax()
        ax.set_ylim(0, max_y)
        ax.set_xlim(0, max_x)
        ax.plot(max_s)
        # ax.set_title(title)
        for path in paths_s:
            plt.savefig(path.format('stress_curve'))
        if show_plot: plt.show()
        else : plt.close()

def figs_stress_scatter(stresses, paths, ordered=True, show_plot=False):
    max_y = 1.05*np.array(stresses).max()
    steps = np.array(stresses).shape[1] -1
    for step, s_ in enumerate(np.array(stresses).transpose((1,0,2))):
        paths_s = [path.format(step/steps*100) for path in paths]
        ax = get_ax()
        ax.set_ylim(0, max_y)
        for s in s_:
            if ordered: s = np.flip(np.sort(s))
            ax.scatter(range(len(s)), s)
        # ax.set_title(title)
        for path in paths_s:
            if ordered: plt.savefig(path.format('ordered_stress_scatter'))
            else: plt.savefig(path.format('stress_scatter'))
        if show_plot: plt.show()
        else : plt.close()

def figs_energy_curve(energies, paths, show_plot=False):
    max_y = 1.05*np.array(energies).max()
    max_x = np.array(energies).shape[1]
    steps = max_x-1
    e = []
    for step, e_ in enumerate(np.array(energies).transpose()):
        paths_s = [path.format(step/steps*100) for path in paths]
        e.append(e_)
        ax = get_ax()
        ax.set_ylim(0, max_y)
        ax.set_xlim(0, max_x)
        ax.plot(e)
        # ax.set_title(title)
        for path in paths_s:
            plt.savefig(path.format('el_energy'))
        if show_plot: plt.show()
        else : plt.close()
    
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
        
def _ax_arms_as_stress(ax, connectivity, s_matrix, min_, max_, positions):
    arms_pos = _get_arms_pos(connectivity, positions)
    colors = _get_arms_color(connectivity, s_matrix, min_, max_)
    for arm,c in zip(arms_pos, colors):
        ax.plot(arm[:,0], arm[:,1], c=c)#, linewidth=8)

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