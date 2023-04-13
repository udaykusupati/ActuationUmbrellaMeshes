import numpy as np

import csv
import sys
import os
sys.path.append('..')
from configuration import deploy_umbrella_pin_rigid_motion

import helpers_tools as help_
import helpers_grid as help_grid

# ======================================================================
# ============================================================ GENERAL =
# ======================================================================

def create_dir_hierarchy(category_name, degree, rows, cols, deployment, folder_name):
    path_name = 'outputs/'+category_name
    if not os.path.exists(path_name): os.makedirs(path_name)
    folder_name = f'{degree:0>2}_{rows:0>2}_{cols:0>2}_{deployment}_{folder_name}'
    path_name += f'/{folder_name}'
    
    # folder to save gif, jpg, png and csv
    help_.create_dir(path_name)
    help_.create_dir(path_name+'/png/gif')
    help_.create_dir(path_name+'/jpg/gif')
    help_.create_dir(path_name+'/results')
    
    return folder_name, path_name

def write_metadata(path, degree, rows, cols, deployment, steps, active_cells, target_percents):
    with open(path+"/metadata.txt", "w") as f:
        f.write("Degree: " + str(degree) + '\n')
        f.write("Rows  : " + str(rows)   + '\n')
        f.write("Cols  : " + str(cols)   + '\n')
        f.write("Deployment : " + str(deployment) + '\n')
        f.write("Steps      : " + str(steps)      + '\n')
        f.write("Active Cells    : " + str(active_cells)    + '\n')
        f.write("Target Percents : " + str(target_percents) + '\n')

def read_metadata(path):
    with open(path+"/metadata.txt") as f:
        metadata = f.read().splitlines()
    degree = int(metadata[0][8:])
    rows = int(metadata[1][8:])
    cols = int(metadata[2][8:])
    deployment = metadata[3][13:].strip()
    steps = int(metadata[4][13:])
    active_cells = eval(metadata[5][18:])
    target_percents = eval(metadata[6][18:])
    return degree, rows, cols, deployment, steps, active_cells, target_percents

def deploy_in_steps(curr_um, input_data, init_heights, plate_thickness, active_cells, target_percents, path_name,
                    steps=10, stress_type='maxBending', verbose=True, dep='linear'):
        dep_weights = help_grid.set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
        
        # deployent in steps
        stresses_per_steps = np.zeros((steps+1, curr_um.numUmbrellas(), curr_um.numUmbrellas())) # 1st step is deployment 0%
        percents_per_steps = []
        heights            = []
        # positions          = []

        path = path_name+'/results'
        path_stresses = path+'/stresses'
        help_.create_dir(path_stresses)

        # write connectivity
        with open(path+'/connectivity.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(input_data['umbrella_connectivity'])
            
        # write initial position
        with open(path+'/position.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(help_grid.get_center_position(curr_um))
        
        for s in range(steps+1):
            if   dep=='linear': target_percents_step = [p*s/steps for p in target_percents]
            elif dep=='min'   : target_percents_step = [min(p, 100*s/steps) for p in target_percents]
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
                stresses_per_steps[s] = help_.get_max_stress_matrix(curr_um, stress_type)
                heights.append(curr_um.umbrellaHeights)
                # positions.append(help_grid.get_center_position(curr_um))
                # write resuts
                with open(path_stresses+f'/{s:0>2}.csv',"w", newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerows(stresses_per_steps[s])
                    if verbose: print(f'step {s: >2}/{steps} saved.')
                
            else: raise ValueError(f'did not converge at step {s}.')

        # write heights
        with open(path+'/heights.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(np.array(heights)-input_data['thickness'])
        # write percents
        with open(path+'/percents.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(active_cells)
            writer.writerows(percents_per_steps)
        
        return stresses_per_steps, percents_per_steps

def read_results(path_name):
    path = path_name+'/results'
    path_stresses = path+'/stresses'

    connectivity  = [[int(i), int(j)] for [i,j] in _read_csv(path+'/connectivity.csv')]
    init_position = _read_csv(path+'/position.csv')

    heights  = _read_csv(path+'/heights.csv')
    actuation = _read_csv(path+'/percents.csv')
    # format data:
    active_cells = [int(c) for c in actuation[0]] # 1st line is the activated cells indexes[0]
    percents_per_steps = actuation[1:]

    stresses = []
    for i in range(len(heights)):
        stresses.append(_read_csv(path_stresses+f'/{i:0>2}.csv'))

    return connectivity,\
           np.array(init_position),\
           heights,\
           active_cells,\
           percents_per_steps,\
           np.array(stresses)



def get_indexes(degree, rows, cols):
    if degree==3:
        indexes = []
        for r in range(rows):
            for c in range(r*2*cols, (r+1)*2*cols, 2):
                if r%2 == 0:
                    indexes.extend([c, c+1])
                else:
                    indexes.extend([c+1, c])
        return indexes
    
    if degree==4:
        return list(range(rows*cols))

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



# ----------------------------------------------------------------------
# ------------------------------------------------------------ helpers -
# ----------------------------------------------------------------------
def _read_csv(path):
    out = []
    with open(path,"r", newline='') as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            out.append(row)
    return out
