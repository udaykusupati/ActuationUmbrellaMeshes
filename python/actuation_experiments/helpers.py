import math
import numpy as np
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

def deploy_in_steps(curr_um, input_data, init_heights, plate_thickness, active_cells, target_percents,
                    steps=10, stress_type='maxBending', verbose=True):
        dep_weights = set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
        
        # deployent in steps
        stresses_per_steps = np.zeros((steps+1, curr_um.numUmbrellas(), curr_um.numUmbrellas())) # 1st step is deployment 0%
        percents_per_steps = []

        for s in range(steps+1):
            target_percents_step = [p*s/steps for p in target_percents]
            
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
           show_height=False, active_cells=[], target_percents=[]):
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
                 active_cells=[], target_percents=[]):
    center_position = get_center_position(curr_um)
    vertices = input_data['umbrella_connectivity']
    center_xy = np.array(list(zip(center_position[:,0], center_position[:,1])))
    segments = center_xy[vertices]
    
    _, ax = plt.subplots()
    for s in segments:
        ax.plot(s[:,0], s[:,1], c='black')
        _ax_dot_active_cell(ax, active_cells, target_percents, center_position)
    ax.axis('equal')
    plt.show()
    
    
def plot2D_steps(input_data, active_cells, percents_per_steps, init_center_pos, stresses_per_steps,
                 stress_type='maxBending', dir_name='00_test',
                 show_percent=False, show_plot=True):
    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()
    max_x, max_y = steps, stresses_all_nz.max()

    # folder to save images
    dir_name = f'./images/{dir_name}'
    _create_dir(dir_name)

    src = np.array(input_data['umbrella_connectivity'])[:,0]
    dst = np.array(input_data['umbrella_connectivity'])[:,1]
    max_stress_per_arm = stresses_per_steps.transpose()[src,dst].max(axis=1)
        
    for s, (s_matrix, percents) in enumerate(zip(stresses_per_steps,percents_per_steps)):
        min_stress_step = s_matrix[s_matrix != 0].min()
        max_stress_step = s_matrix[s_matrix != 0].max()
        if s==0: min_stress_step, max_stress_step = 0,0 # manage random perturbation, at step 0, no deploymen at all
        max_stresses.append(max_stress_step)

        title = f'{s/steps*100:.0f}% deployed\n{stress_type}: {max_stress_step:.2f}'
        fig_saved_name = f'{dir_name}/{stress_type}_'+'{}'+f'_{s/steps*100:0>3.0f}Deployment.jpg'

        # normalized with general extrems values
        _fig_arm_stresses(input_data, active_cells, percents, init_center_pos, show_percent,
                          s_matrix, min_stress_all, max_stress_all, show_plot, title, fig_saved_name.format('structure_all'))
        
        # normalized with step extrems values
        _fig_arm_stresses(input_data, active_cells, percents, init_center_pos, show_percent,
                          s_matrix, min_stress_step, max_stress_step, show_plot, title, fig_saved_name.format('structure_perSteps'))

        # normalized with own extrems values
        _fig_stress_compare_own(input_data, active_cells, percents, init_center_pos, show_percent, s_matrix, max_stress_per_arm, show_plot, title,
                                fig_saved_name.format('structure_own'))
        
        # stress curve
        _fig_stress_curve(max_stresses, max_x, max_y, show_plot, title, fig_saved_name.format('sPlot'))
        
    
# ------------------------------------------------------------ helpers -
def _ax_dot_active_cell(ax, active_cells, target_percents, positions):
    for i, p in zip(active_cells, target_percents):
        r = p/100
        [x,y,_] = positions[i]
        ax.scatter(x,y, color=(r,1-r,0))

def _ax_arms_as_stress(ax, input_data, s_matrix, min_, max_, position):
    for i,j in input_data['umbrella_connectivity']:
        c = _color_map(s_matrix[i,j], min_, max_)
        [x,y,_] = position[[i,j]].transpose()
        ax.plot(x, y, c=c)

def _ax_show_percent(ax, show_percent, active_cells, target_percents, position): #### HERER : Add ax as param, change in experience function (steps)
    for i, p in zip(active_cells*show_percent, target_percents):
        ax.annotate(f'{p: >5.0f}', position[i][:2], ha='left', va='center', color='black')

def _ax_plot_stresses(ax, input_data, stress_matrix, min_, max_, active_cells, percents, position, show_percent):    
    _ax_arms_as_stress(ax, input_data, stress_matrix, min_, max_, position)
    _ax_dot_active_cell(ax, active_cells, percents, position)
    _ax_show_percent(ax, show_percent, active_cells, percents, position)
    
def _fig_arm_stresses(input_data, active_cells, percents, init_center_pos, show_percent, s_matrix, min_, max_, show_plot, title, file_name=''):
    fig_size = 8
    _, ax = plt.subplots(figsize=(fig_size, fig_size))
    _ax_plot_stresses(ax, input_data, s_matrix, min_, max_, active_cells, percents, init_center_pos, show_percent)
    ax.set_title(title)
    ax.axis('equal')
    if file_name != '': plt.savefig(file_name)
    if show_plot: plt.show()
    plt.close()
    
def _fig_stress_compare_own(input_data, active_cells, percents, position, show_percent, s_matrix, max_ls, show_plot, title, file_name=''):
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
    if file_name != '': plt.savefig(file_name)
    if show_plot: plt.show()
    plt.close()

def _fig_stress_curve(max_stresses, max_x, max_y, show_plot, title, file_name=''):
    fig_size = 8
    _, ax = plt.subplots(figsize=(fig_size, fig_size))
    ax.plot(max_stresses)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)
    ax.set_title(title)
    if file_name != '': plt.savefig(file_name)
    if show_plot: plt.show()
    plt.close()

def _color_map(value, min_, max_, expo=1):
    if max_ == min_:
        if max_ == 0: return 0,1,0 # no stress at all -> green
        else:         return 0,0,0 # max stress/height everywhere
    r = ((value-min_)/(max_-min_))**expo # [xRK] thy some expo <1
    g = 1-r
    b = 0
    return r,g,b

def _get_smatrix_min_max(curr_um, stress_type,
                         zero_as_extrem=False):
    matrix = _get_max_stress_matrix(curr_um, stress_type)
    if not zero_as_extrem:
        tmp = matrix[np.nonzero(matrix)]
        return matrix, tmp.min(), tmp.max()
    return matrix, matrix.min(), matrix.max()

def _create_dir(name):
    if not os.path.exists(name): os.makedirs(name)
    else: raise ValueError(f'folder {name} already exists')


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