import numpy as np

import csv
import sys
import os
sys.path.append('..')
from configuration import deploy_umbrella_pin_rigid_motion
from pipeline_helper import allEnergies

import helpers_tools as help_
import helpers_grid as help_grid
from plots import projection2D

# ======================================================================
# ============================================================ GENERAL =
# ======================================================================

def create_dir_hierarchy(category_name, degree, rows, cols, deployment, folder_name):
    folder_name = f'{degree:0>2}_{rows:0>2}_{cols:0>2}_{folder_name}'
    path = 'outputs/'+category_name + f'/{folder_name}'
    if not os.path.exists(path): os.makedirs(path)

    path_dep = path+f'/{deployment}_deployment'
    if not os.path.exists(path_dep): os.makedirs(path_dep)
    else: raise ValueError(f'deployment {deployment} already computed.')

    for type_ in ['heights', 'energies']:
        _sub_folder(path_dep+f'/{type_}')

    for s_type in help_.get_stresses_types():
        _sub_folder(path_dep+f'/stresses/{s_type}')
    
    return folder_name, path

def write_metadata(path, degree, rows, cols, deployment, steps, active_cells, target_percents):
    with open(path+"/metadata.txt", "w") as f:
        f.write("Degree: " + str(degree) + '\n')
        f.write("Rows  : " + str(rows)   + '\n')
        f.write("Cols  : " + str(cols)   + '\n')
        f.write("Deployment : " + str(deployment) + '\n')
        f.write("Steps      : " + str(steps)      + '\n')
        f.write("Active Cells    : " + str(active_cells)    + '\n')
        f.write("Target Percents : " + str(target_percents) + '\n')

def deploy_in_steps(curr_um, input_data, init_heights, plate_thickness, active_cells, target_percents, path, deployment,
                    steps=10, verbose=True, dep='linear'):
        dep_weights = help_grid.set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
        
        # deployent in steps
        stresses_per_steps = np.zeros((steps+1, curr_um.numUmbrellas(), curr_um.numUmbrellas())) # 1st step is deployment 0%
        percents_per_steps = []
        heights            = []
        energies           = []

        # header for energies.csv
        energies.append(allEnergies(curr_um).keys()) # will be read back as dict

        # write connectivity
        with open(path+'/connectivity.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(input_data['umbrella_connectivity'])
            
        # write initial position
        with open(path+'/position.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(help_grid.get_center_position(curr_um))
        
        stresses_types = help_.get_stresses_types()
        for s in range(steps+1):
            if   dep=='linear'      : target_percents_step = [p*s/steps for p in target_percents]
            elif dep=='incremental' : target_percents_step = [min(p, 100*s/steps) for p in target_percents]
            elif dep=='max':
                max_p = max(target_percents)
                target_percents_step = [min(p, max_p*s/steps) for p in target_percents]
            else : raise ValueError(f'deployment unknown: {dep} (choises:\'linear\',\'min\',\'max\')')
            
            percents_per_steps.append(target_percents_step)
            target_heights = help_grid.percent_to_height(init_heights, plate_thickness, active_cells, target_percents_step)
            target_height_multiplier = help_grid.set_target_height(curr_um.numUmbrellas(), active_cells, target_heights)
            success, _ = deploy_umbrella_pin_rigid_motion(curr_um,
                                                          plate_thickness,
                                                          target_height_multiplier,
                                                          dep_weights=dep_weights)
            if success:
                heights.append(curr_um.umbrellaHeights)
                energies.append(list(allEnergies(curr_um).values()))
                for s_type in stresses_types:
                    path_stresses = path+f'/{dep}_deployment/stresses/{s_type}/values'
                    stresses_per_steps[s] = help_.get_max_stress_matrix(curr_um, s_type)
                    # write resuts
                    with open(path_stresses+f'/step_{s:0>2}.csv',"w", newline='') as csvfile:
                        writer = csv.writer(csvfile)
                        writer.writerows(stresses_per_steps[s])
                if verbose: print(f'step {s: >2}/{steps} saved.')
                
            else: raise ValueError(f'did not converge at step {s}.')
        
        projection2D(input_data['umbrella_connectivity'], curr_um, active_cells, target_percents, file_name=path+'/projection2D.png')
        
        # write heights
        with open(path+f'/{dep}_deployment/heights/values/heights.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(np.array(heights)-input_data['thickness'])
        # write percents and active cells
        with open(path+f'/{dep}_deployment/heights/values/percents.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(active_cells)
            writer.writerows(percents_per_steps)
        # write energies
        with open(path+f'/{dep}_deployment/energies/values/energies.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(energies)

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


def _sub_folder(path):
    os.makedirs(path+'/values')
    os.makedirs(path+'/png/gif')
    os.makedirs(path+'/jpg/gif')
