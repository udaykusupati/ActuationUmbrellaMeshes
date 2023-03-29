import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('..')
import umbrella_mesh
import linkage_vis
from visualization_helper import get_color_field

## ===============
## PLOT 2D - close Umbrella
## ===============

def get_center_position(curr_um):
    nb_cell = curr_um.numUmbrellas()
    center_position = np.zeros([nb_cell, 3])
    for i in range(nb_cell):
        top_idx = curr_um.getUmbrellaCenterJi(i, 0)
        center_position[i] = curr_um.joint(top_idx).position
    return center_position

def plot2D(input_data, curr_um, show_height=False, active_cells=[], target_percents=[]):
    center_position = get_center_position(curr_um)
    vertices = input_data['base_mesh_v']
    edge     = np.roll(np.insert(input_data['base_mesh_f'], 0, input_data['base_mesh_f'][:,0], axis=1), -1, axis=1)

    ## == plotting == ##
    fig_length = 15#max(5, len(center_position)**0.5)
    fig = plt.figure(figsize=(fig_length, fig_length))
    # cells' edges
    for e in vertices[edge]:
        plt.plot(e[:, 0], e[:, 1], color="lightblue")
    # cell index at its center
    for i, [x,y,z] in enumerate(center_position):
        plt.annotate(f'{i}', (x,y), ha='center', color='black')
    if show_height:
        h_max, h_min = max(curr_um.umbrellaHeights), min(curr_um.umbrellaHeights)
        for [x,y,z], h in zip(center_position, curr_um.umbrellaHeights):
            if h_max != h_min :
                r = (h-h_min)/(h_max-h_min)
                g = 0.1
                b = 1-r
            else:
                r,g,b = 0,0,0
            plt.annotate(f'[{h:.2f}]', (x,y), textcoords='offset points', xytext=(0,-12), ha='center', c=(r,g,b))
    # overwrite index if actve_cell
    if active_cells != []:
        for i, p in zip(active_cells, target_percents):
            g = p/100
            r = 1-g
            b = 0
            [x,y,_] = center_position[i]
            plt.annotate(f'{i}', (x,y), ha='center', color=(r,g,b), weight='bold')
    
    # show plot
    plt.axis('equal')
    plt.show()
    
def projection2D(input_data, curr_um, active_cells=[], target_percents=[]):
    center_position = get_center_position(curr_um)
    vertices = input_data['umbrella_connectivity']
    center_xy = np.array(list(zip(center_position[:,0], center_position[:,1])))
    segments = center_xy[vertices]
    
    for s in segments:
        plt.plot(s[:,0], s[:,1], c='black')
    
    # color for actve_cell
    if active_cells != []:
        for i, p in zip(active_cells, target_percents):
            r = p/100
            g = 1-r
            b = 0
            [x,y,_] = center_position[i]
            plt.scatter(x,y, color=(r,g,b))
            
    plt.axis('equal')
    plt.show()


## ==============
## PLOT 3D
## ===============
def plot3D(curr_um, input_data, rod_colors=None, uidBased=False):
    if rod_colors is None:
        rod_colors = get_color_field(curr_um, input_data, uidBased) 
    view = linkage_vis.LinkageViewer(curr_um, width=800, height=600)
    view.update(scalarField = rod_colors)
    return view
    
## ===============
## Constraints
## ===============
def set_actives_dep_weights(numUmbrellas, *active_cells, dep_factors = np.logspace(-4, 0, 5)):
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

def plot_stress3D(curr_um, stress_type, verbose=True):
    view = linkage_vis.LinkageViewer(curr_um, width=800, height=600)
    stresses = _set_plate_stress_null(curr_um, _get_stresses(curr_um, stress_type))
    view.update(scalarField = stresses)
    if verbose:
        print(f'{stress_type} Stresses Extrem values:\n\
    max : {np.array(stresses).max():.2e}\n\
    min : {np.array(stresses).min():.2e}')
    return view

def _get_stresses(curr_um, stress_type):
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

def _set_plate_stress_null(curr_um, stresses):
    stresses = np.array(stresses)
    for sid in range(curr_um.numSegments()):
        seg = curr_um.segment(sid)
        if seg.segmentType() == umbrella_mesh.SegmentType.Plate:
            stresses[sid] = 0
    return stresses.tolist()


def _get_stress_matrix(grid, stress_type='maxBending'):
    '''return the adjacency max stress matrix'''
    matrix = np.zeros((grid.numUmbrellas, grid.numUmbrellas))
    stress_all = _get_stresses(grid.curr_um, stress_type)
    for seg_id in range(grid.curr_um.numSegments()):
        seg = grid.curr_um.segment(seg_id)
        if seg.segmentType()==umbrella_mesh.SegmentType.Arm:
            sjoint = seg.startJoint
            ejoint = seg.endJoint
            stresses = stress_all[seg_id]
            if grid.curr_um.joint(sjoint).jointPosType()==umbrella_mesh.JointPosType.Arm:
                i,j = grid.curr_um.joint(ejoint).umbrellaID()
            if grid.curr_um.joint(ejoint).jointPosType()==umbrella_mesh.JointPosType.Arm:
                i,j = grid.curr_um.joint(sjoint).umbrellaID()
            matrix[i,j] = matrix[j,i] = max(matrix[i,j], max(stresses))
    return matrix

def plot_stress2D(grid, active_cells=[], target_percents=[], stress_type='maxBending', zero_as_extrem = False, show_percent=False):
    matrix = _get_stress_matrix(grid, stress_type)
    max_, min_ = matrix.max(), matrix.min()
    if not zero_as_extrem:
        matrix_nz = matrix[np.nonzero(matrix)]
        max_, min_ = matrix_nz.max(),matrix_nz.min()
    color = []
    for i,j in grid.input_data['umbrella_connectivity']:
        if max_ != min_:
            r = (matrix[i,j]-min_)/(max_-min_)
            g = 1-r
            b = 0
        else: r,g,b = 0,0,0
        [x,y,_] = grid.init_center_xyz[[i,j]].transpose()
        plt.plot(x, y, c=(r,g,b))
    
    # color for actve_cell
    if active_cells != []:
        for i, p in zip(active_cells, target_percents):
            r = p/100
            g = 1-r
            b = 0
            [x,y,_] = grid.init_center_xyz[i]
            plt.scatter(x,y, color=(r,g,b))
            if show_percent:
                plt.annotate(f'{p: >5}',(x,y), ha='left', va='center', color='black')
            
    plt.axis('equal')
    plt.show()
    
'''  
def projection2D(input_data, curr_um, active_cells=[], target_percents=[]):
    center_position = get_center_position(curr_um)
    vertices = input_data['umbrella_connectivity']
    center_xy = np.array(list(zip(center_position[:,0], center_position[:,1])))
    segments = center_xy[vertices]
    
    for s in segments:
        plt.plot(s[:,0], s[:,1], c='black')
    
    # color for actve_cell
    if active_cells != []:
        for i, p in zip(active_cells, target_percents):
            g = p/100
            r = 1-g
            b = 0
            [x,y,_] = center_position[i]
            plt.scatter(x,y, color=(r,g,b))
            
    plt.axis('equal')
    plt.show()
 '''