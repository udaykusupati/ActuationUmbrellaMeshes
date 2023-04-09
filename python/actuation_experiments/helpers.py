import math
import numpy as np
import matplotlib as mpl
from matplotlib import colormaps
import matplotlib.pyplot as plt

import sys
import os
sys.path.append('..')
import umbrella_mesh
import linkage_vis
from visualization_helper import get_color_field
from configuration import deploy_umbrella_pin_rigid_motion

# ======================================================================
# ============================================================ GENERAL =
# ======================================================================

def create_dir_hierarchy(category_name, degree, rows, cols, deployment, folder_name):
    path_name = 'outputs/'+category_name
    if not os.path.exists(path_name): os.makedirs(path_name)
    folder_name = f'{degree:0>2}_{rows:0>2}_{cols:0>2}_{deployment}_{folder_name}'
    path_name += f'/{folder_name}'
    
    # folder to save gif, jpg, png and csv
    _create_dir(path_name)
    _create_dir(path_name+'/png/gif')
    _create_dir(path_name+'/jpg/gif')
    _create_dir(path_name+'/csv')
    
    return folder_name, path_name


def _create_dir(name):
    if not os.path.exists(name): os.makedirs(name)
    else: raise ValueError(f'folder {name} already exists')


def get_center_position(curr_um):
    nb_cell = curr_um.numUmbrellas()
    center_position = np.zeros([nb_cell, 3])
    for i in range(nb_cell):
        top_idx = curr_um.getUmbrellaCenterJi(i, 0)
        center_position[i] = curr_um.joint(top_idx).position
    return center_position

def set_actives_dep_weights(numUmbrellas, *active_cells,
                            dep_factors = np.logspace(-4, 0, 5)):
    weights_per_cell = np.zeros(numUmbrellas)
    weights_per_cell[active_cells] = 1
    return np.einsum('i, j -> ij', dep_factors, weights_per_cell)

def set_target_height(numUmbrellas, active_cell, target_heights):
    target_height_multiplier = np.ones(numUmbrellas)*-1
    target_height_multiplier[active_cell] = target_heights
    return target_height_multiplier

def percent_to_height(init_height, thickness, indexes, percents):
    return [((1-percent/100)*(init_height[idx]-thickness)+thickness)/thickness
            for percent,idx in zip(percents,indexes)]

def linear_heights(a,b, step=1):
    '''
    from undeployed to deployed :
    (0,4) -> [0,1,2,3,4],[100, 75, 50, 25, 0]
    (4,0) -> [0,1,2,3,4],[0, 25, 50, 75, 100]
    '''
    active_cells    = []
    target_percents = []
    if b>a:
        len_ = b-a
        for i in range(0, step*len_+1,step):
            active_cells.append(a+i)
            target_percents.append(100*i/len_)
    elif a>b:
        len_ = a-b
        for i in range(0, step*len_+1, step):
            active_cells.append(step*len_+b-i)
            target_percents.append(100*i/len_)

    return active_cells, target_percents

# def linear_ 
### --> DECOMPOSE FUNCTION -> give lis to receive height (for diagonal)
 
def linear_height_ls_idx(ls, min_dep=0, max_dep=100):
    max_ = max(ls)
    min_ = min(ls)
    percents = []
    for idx in ls:
        p = (idx-min_)/(max_-min_)
        percents.append(min_dep+p*max_dep)
    return percents

def linear_height_ls(ls, min_dep=0, max_dep=100):
    len_ = len(ls)
    percents = []
    for idx in range(len_):
        p = idx/(len_-1)
        percents.append(min_dep+p*max_dep)
    return percents

def deploy_in_steps(curr_um, input_data, init_heights, plate_thickness, active_cells, target_percents,
                    steps=10, stress_type='maxBending', verbose=True, dep='linear'):
        dep_weights = set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
        
        # deployent in steps
        stresses_per_steps = np.zeros((steps+1, curr_um.numUmbrellas(), curr_um.numUmbrellas())) # 1st step is deployment 0%
        percents_per_steps = []

        for s in range(steps+1):
            if   dep=='linear': target_percents_step = [p*s/steps for p in target_percents]
            elif dep=='min'   : target_percents_step = [min(p, 100*s/steps) for p in target_percents]
            elif dep=='max':
                max_p = max(target_percents)
                target_percents_step = [min(p, max_p*s/steps) for p in target_percents]
            else : raise ValueError(f'deployment unknown: {dep} (choises:\'linear\',\'min\',\'max\')')
            
            percents_per_steps.append(target_percents_step)
            target_heights = percent_to_height(init_heights, plate_thickness, active_cells, target_percents_step)
            target_height_multiplier = set_target_height(curr_um.numUmbrellas(), active_cells, target_heights)
            success, _ = deploy_umbrella_pin_rigid_motion(curr_um,
                                                          plate_thickness,
                                                          target_height_multiplier,
                                                          dep_weights=dep_weights)
            if success:
                stresses_per_steps[s] = _get_max_stress_matrix(curr_um, stress_type)
                if verbose: print(f'step {s: >2}/{steps} computed.')
            else: raise ValueError(f'did not converge at step {s}.')
        
        return stresses_per_steps, percents_per_steps


# ======================================================================
# ============================================================== PLOTS =
# ======================================================================
# ------------------- 3D -------------------
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
    stresses = _set_plate_stress_null(curr_um, _get_stresses(curr_um, stress_type))
    view.update(scalarField = stresses)
    if verbose:
        print(f'{stress_type} Stresses Extrem values:\n\
    max : {np.array(stresses).max():.2e}\n\
    min : {np.array(stresses).min():.2e}')
    return view

# ------------------- 2D -------------------
def plot2D(input_data, curr_um,
           show_height=False, active_cells=[], target_percents=[], file_name=''):
    center_position = get_center_position(curr_um)
    vertices = input_data['base_mesh_v']
    edge     = np.roll(np.insert(input_data['base_mesh_f'], 0,
                                 input_data['base_mesh_f'][:,0], axis=1),
                       -1, axis=1)
    # plot:
    fig_size = 15
    _, ax = plt.subplots(figsize=(fig_size, fig_size))
    # edges
    for e in vertices[edge]:
        ax.plot(e[:, 0], e[:, 1], color="lightblue")
    # index
    for i, [x,y,z] in enumerate(center_position):
        ax.annotate(f'{i}', (x,y), ha='center', color='black')
    # height
    um_heights = curr_um.umbrellaHeights
    for [x,y,z], h in zip(show_height*center_position, um_heights):
        c = _color_map(h, min(um_heights), max(um_heights))
        ax.annotate(f'[{h:.2f}]', (x,y), textcoords='offset points', xytext=(0,-12), ha='center', c=c)
    # actve cells
    for i, p in zip(active_cells, target_percents):
        r = p/100
        ax.annotate(f'{i}', center_position[i][:2], ha='center', color=(r,1-r,0), weight='bold')
    # show plot
    ax.axis('equal')
    if file_name != '':
        plt.savefig(file_name)
    plt.show()

def plot2D_stress(curr_um, input_data, init_center_pos,
                  active_cells=[], target_percents=[], stress_type='maxBending',
                  zero_as_extrem = False, show_percent=False):
    stress_matrix, min_, max_ = _get_smatrix_min_max(curr_um, stress_type, zero_as_extrem)
    
    _, ax = plt.subplots()
    _ax_plot_stresses(ax, input_data, stress_matrix, min_, max_, active_cells, target_percents, init_center_pos, show_percent)
    ax.axis('equal')
    plt.show()
    
def projection2D(input_data, curr_um,
                 active_cells=[], target_percents=[], file_name=''):
    center_position = get_center_position(curr_um)
    vertices = input_data['umbrella_connectivity']
    center_xy = np.array(list(zip(center_position[:,0], center_position[:,1])))
    segments = center_xy[vertices]
    
    _, ax = plt.subplots()
    for s in segments:
        ax.plot(s[:,0], s[:,1], c='black')
        _ax_dot_active_cell(ax, active_cells, target_percents, center_position)
    ax.axis('equal')
    if file_name != '':
        plt.savefig(file_name)
    plt.show()
    
    
def plot2D_steps(input_data, active_cells, percents_per_steps, init_center_pos, stresses_per_steps,
                 stress_type='maxBending', dir_name='00_test',
                 show_percent=False, show_plot=True):
    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()
    max_x, max_y = steps, stresses_all_nz.max()+5

    src = np.array(input_data['umbrella_connectivity'])[:,0]
    dst = np.array(input_data['umbrella_connectivity'])[:,1]
    max_stress_per_arm = stresses_per_steps.transpose()[src,dst].max(axis=1)
    
    title = '{:.0f}% deployed\n'+f'{stress_type}:'+' {:.2f}'
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
        _fig_arm_stresses(input_data, active_cells, percents, init_center_pos, show_percent,
                          s_matrix, deployed*min_stress_all, deployed*max_stress_all, show_plot, title_s, path_names_s, 'structure_all')
        
        # normalized with step extrems values
        _fig_arm_stresses(input_data, active_cells, percents, init_center_pos, show_percent,
                          s_matrix, deployed*min_stress_step, deployed*max_stress_step, show_plot, title_s, path_names_s, 'structure_perSteps')

        # normalized with own extrems values
        _fig_stress_compare_own(input_data, active_cells, percents, init_center_pos, show_percent, s_matrix, deployed*max_stress_per_arm, show_plot,
                                title_s, path_names_s, 'structure_own')
        
        # ordered stresses
        _fig_stress_scatter(deployed*s_matrix, max_y, show_plot, title_s, path_names_s, 'scatter')
        
        # stress curve
        _fig_stress_curve(deployed*max_stresses, max_x, max_y, show_plot, title_s, path_names_s, 'sPlot')
        
        # perturbations do not affect deployed state
        deployed = True
        
    
# ------------------------------------------------------------ helpers -
def _ax_dot_active_cell(ax, active_cells, target_percents, positions):
    for i, p in zip(active_cells, target_percents):
        r = p/100
        [x,y,_] = positions[i]
        ax.scatter(x,y, color=(r,1-r,0))

def _ax_arms_as_stress(ax, input_data, s_matrix, min_, max_, position):
    for i,j in input_data['umbrella_connectivity']:
        c = _color_map(s_matrix[i,j], min_, max_, '') # turbo
        [x,y,_] = position[[i,j]].transpose()
        ax.plot(x, y, c=c)

def _ax_show_percent(ax, show_percent, active_cells, target_percents, position): #### HERER : Add ax as param, change in experience function (steps)
    for i, p in zip(active_cells*show_percent, target_percents):
        ax.annotate(f'{p: >5.0f}', position[i][:2], ha='left', va='center', color='black')

def _ax_plot_stresses(ax, input_data, stress_matrix, min_, max_, active_cells, percents, position, show_percent):    
    _ax_arms_as_stress(ax, input_data, stress_matrix, min_, max_, position)
    _ax_dot_active_cell(ax, active_cells, percents, position)
    _ax_show_percent(ax, show_percent, active_cells, percents, position)
    
def _fig_arm_stresses(input_data, active_cells, percents, init_center_pos, show_percent, s_matrix, min_, max_, show_plot, title, path_names, file_name=''):
    fig_size = 8
    _, ax = plt.subplots(figsize=(fig_size, fig_size))
    _ax_plot_stresses(ax, input_data, s_matrix, min_, max_, active_cells, percents, init_center_pos, show_percent)
    ax.set_title(title)
    ax.axis('equal')
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()
    
def _fig_stress_compare_own(input_data, active_cells, percents, position, show_percent, s_matrix, max_ls, show_plot, title, path_names, file_name=''):
    fig_size = 8
    _, ax = plt.subplots(figsize=(fig_size, fig_size))

    for arm,(i,j) in enumerate(input_data['umbrella_connectivity']):
        c = _color_map(s_matrix[i,j], 0, max_ls[arm]) # min is always 0s (undeployed state)
        [x,y,_] = position[[i,j]].transpose()
        ax.plot(x, y, c=c)
    _ax_dot_active_cell(ax, active_cells, percents, position)
    _ax_show_percent(ax, show_percent, active_cells, percents, position)
    
    ax.set_title(title)
    ax.axis('equal')
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()

def _fig_stress_curve(max_stresses, max_x, max_y, show_plot, title, path_names, file_name=''):
    fig_size = 8
    _, ax = plt.subplots(figsize=(fig_size, fig_size))
    ax.plot(max_stresses)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)
    ax.set_title(title)
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()

def _fig_stress_scatter(s_matrix, max_y, show_plot, title, path_names, file_name=''):
    fig_size = 8
    _, ax = plt.subplots(figsize=(fig_size, fig_size))
    sort = np.flip(np.unique(s_matrix[s_matrix!=0]))
    ax.scatter(range(len(sort)), sort, s=5)
    ax.set_ylim(0, max_y)
    ax.set_title(title)
    for path in path_names:
        plt.savefig(path.format(file_name))
    if show_plot: plt.show()
    plt.close()
    
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

def _get_smatrix_min_max(curr_um, stress_type,
                         zero_as_extrem=False):
    matrix = _get_max_stress_matrix(curr_um, stress_type)
    if not zero_as_extrem:
        tmp = matrix[np.nonzero(matrix)]
        return matrix, tmp.min(), tmp.max()
    return matrix, matrix.min(), matrix.max()

# ======================================================================
# ============================================================ HELPERS =
# ======================================================================
def _set_plate_stress_null(curr_um,
                           stresses):
    stresses = np.array(stresses)
    for sid in range(curr_um.numSegments()):
        seg = curr_um.segment(sid)
        if seg.segmentType() == umbrella_mesh.SegmentType.Plate:
            stresses[sid] = 0
    return stresses.tolist()

def _get_max_stress_matrix(curr_um,
                       stress_type='maxBending'):
    '''return the adjacency max stress matrix'''
    matrix = np.zeros((curr_um.numUmbrellas(), curr_um.numUmbrellas()))
    stress_all = _get_stresses(curr_um, stress_type)
    for seg_id in range(curr_um.numSegments()):
        seg = curr_um.segment(seg_id)
        if seg.segmentType()==umbrella_mesh.SegmentType.Arm:
            sjoint = seg.startJoint
            ejoint = seg.endJoint
            stresses = stress_all[seg_id]
            if curr_um.joint(sjoint).jointPosType()==umbrella_mesh.JointPosType.Arm:
                i,j = curr_um.joint(ejoint).umbrellaID()
            if curr_um.joint(ejoint).jointPosType()==umbrella_mesh.JointPosType.Arm:
                i,j = curr_um.joint(sjoint).umbrellaID()
            matrix[i,j] = matrix[j,i] = max(matrix[i,j], max(stresses))
    return matrix

def _get_stresses(curr_um,
                  stress_type):
    if stress_type=='Von Mises':
        stresses = curr_um.maxVonMisesStresses()
    elif stress_type=='maxBending':
        stresses = curr_um.maxBendingStresses()
    elif stress_type=='minBending':
        stresses = curr_um.minBendingStresses()
    elif stress_type=='Twisting':
        stresses = curr_um.twistingStresses()
    elif stress_type=='Stretching':
        stresses = curr_um.stretchingStresses()
    else:
        raise ValueError(f'the required stress type <{stress_type}> do not correspond to any available stress')
    return stresses
            