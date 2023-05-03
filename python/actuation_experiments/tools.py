import numpy as np
import glob
from PIL import Image

import csv
import sys
import os
sys.path.append('..')
from configuration import deploy_umbrella_pin_rigid_motion
from pipeline_helper import allEnergies

import helpers_tools as help_
import helpers_grid as help_grid

# ======================================================================
# =============================================================== GRID =
# ======================================================================
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

def get_center_position(curr_um):
    nb_cell = curr_um.numUmbrellas()
    center_position = np.zeros([nb_cell, 3])
    for i in range(nb_cell):
        top_idx = curr_um.getUmbrellaCenterJi(i, 0)
        center_position[i] = curr_um.joint(top_idx).position
    return center_position

# ======================================================================
# ============================================================ GENERAL =
# ======================================================================
def create_dir_hierarchy(category_name, degree, rows, cols, deployment, folder_name):
    folder_name = f'{degree:0>2}_{rows:0>2}_{cols:0>2}_{folder_name}' if rows!=0 and cols!=0 \
             else f'{degree:0>2}_{folder_name}'
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

def write_metadata(path, folder_name, degree, rows, cols, steps, active_cells_per_phase, target_percents_per_phase):
    with open(path+"/metadata.txt", "w") as f:
        f.write("Name  : " + folder_name + '\n')
        f.write("Degree: " + str(degree) + '\n')
        f.write("Rows  : " + str(rows)   + '\n')
        f.write("Cols  : " + str(cols)   + '\n')
        f.write("Steps : " + str(steps)      + '\n')
        f.write("Phases: " + str(len(active_cells_per_phase)) + '\n')
        f.write("Active Cells   : " + str(active_cells_per_phase)    + '\n')
        f.write("Target Percents: " + str(target_percents_per_phase) + '\n')
    

def deploy_in_steps(curr_um, input_data, init_heights, plate_thickness, active_cells_all, target_percents_all, deployment, path,
                    steps=10, verbose=True):

    prev_target_percent = [0]*len(target_percents_all[0])
    max_phase = len(active_cells_all)-1
    for phase, (active_cells, target_percents) in enumerate(zip(active_cells_all, target_percents_all)):
        if verbose:
            print(f'==== PHASE {phase+1: >2} ====')
            print(f'{active_cells    = }')
            print(f'{target_percents = }')
        dep_weights = help_grid.set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
        # deployent in steps
        stresses_per_steps = np.zeros((steps+1, curr_um.numUmbrellas(), curr_um.numUmbrellas())) # 1st step is deployment 0%
        percents_per_steps = []
        heights            = []
        energies           = []

        # header for energies.csv
        energies.append(allEnergies(curr_um).keys()) # will be read back as dict
        
        # write 
        _write_rows(path+'/connectivity.csv',input_data['umbrella_connectivity'])
        _write_rows(path+'/position.csv', get_center_position(curr_um))

        stresses_types = help_.get_stresses_types()
        for s in range(steps+1):
            if   deployment=='linear'      : target_percents_step = [max(0, min(prev_p+(p-prev_p)*s/steps, 100)) for prev_p, p in zip(prev_target_percent, target_percents)]
            elif deployment=='incremental' : target_percents_step = [min(p, 100*s/steps) for p in target_percents]
            elif deployment=='max':
                max_p = max(target_percents)
                target_percents_step = [min(p, max_p*s/steps) for p in target_percents]
            else : raise ValueError(f'deployment unknown: {deployment} (choises:\'linear\',\'incremental\',\'max\')')
            
            percents_per_steps.append(target_percents_step)
            target_heights = help_grid.percent_to_height(init_heights, plate_thickness, active_cells, target_percents_step)
            target_height_multiplier = help_grid.set_target_height(curr_um.numUmbrellas(), plate_thickness, active_cells, target_heights)
            success, _ = deploy_umbrella_pin_rigid_motion(curr_um,
                                                            plate_thickness,
                                                            target_height_multiplier,
                                                            dep_weights=dep_weights)
            if success:
                heights.append(curr_um.umbrellaHeights)
                energies.append(list(allEnergies(curr_um).values()))
                for s_type in stresses_types:
                    path_stresses = path+f'/{deployment}_deployment/stresses/{s_type}/values'
                    stresses_per_steps[s] = help_.get_max_stress_matrix(curr_um, s_type)
                    # write resuts
                    with open(path_stresses+f'/phase{phase+1}_step{s:0>2}.csv',"w", newline='') as csvfile:
                        writer = csv.writer(csvfile)
                        writer.writerows(stresses_per_steps[s])
                if verbose: print(f'step {s: >2}/{steps} saved.')

            else: raise ValueError(f'did not converge at step {s}.')

        # write
        _write_rows(path+f'/{deployment}_deployment/heights/values/phase{phase+1}_heights.csv',
                    np.array(heights)-input_data['thickness'])
        _write_rows(path+f'/{deployment}_deployment/energies/values/phase{phase+1}_energies.csv',
                    energies)
        # write percents and active cells
        with open(path+f'/{deployment}_deployment/heights/values/phase{phase+1}_percents.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(active_cells)
            writer.writerows(percents_per_steps)
        
        if phase < max_phase:
            prev_target_percent = help_grid.height_to_percent(init_heights, plate_thickness, active_cells_all[phase+1], curr_um.umbrellaHeights)

def img_to_gif(path, deployment, stress_type, duration=500, loop=2, verbose=False):
    path = f'{path}/{deployment}_deployment'
    img_to_gif_1D(path, stress_type, duration=duration, loop=loop)
    if verbose: print('gif 1D done')
    img_to_gif_2D(path, stress_type, duration=duration, loop=loop)
    if verbose: print('gif 2D done')

def img_to_gif_1D(path, stress_type, duration=500, loop=2):
    _img_to_gif_stress_1D (path, stress_type, duration=duration, loop=loop)
    _img_to_gif_heights_1D(path,              duration=duration, loop=loop)
    _img_to_gif_energies  (path,              duration=duration, loop=loop)
    
def img_to_gif_2D(path, stress_type, duration=500, loop=2):
    _img_to_gif_stress_2D (path, stress_type, duration=duration, loop=loop)
    _img_to_gif_heights_2D(path,              duration=duration, loop=loop)
    
def img_to_gif_all(path, gif_name, duration=500, loop=2):
    '''
    create gif `gif_name` with all files at path
    '''
    _create_gif(f'{path}/*.*', f'{path}/{gif_name}.gif', duration, loop)

        
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


def _create_gif(images_name, gif_name, duration, loop):
    frames = [Image.open(image) for image in sorted(glob.glob(images_name))]
    frame_one = frames[0]
    frame_one.save(gif_name, append_images=frames[1:], save_all=True, duration=duration, loop=loop)
    
def _gif_list(path, ls, duration, loop):
    for name in ls:
        _create_gif(f'{path}/jpg/phase?_{name}*.jpg',
                   f'{path}/jpg/gif/{name}.gif',
                   duration, loop)
def _img_to_gif_stress_1D(path, stress_type, duration=500, loop=2):
    path_base = f'{path}/stresses/{stress_type}'
    names = ['stress_curve','ordered_stress_scatter', 'stress_scatter']
    _gif_list(path_base, names, duration, loop)
def _img_to_gif_stress_2D(path, stress_type, duration=500, loop=2):
    path_base = f'{path}/stresses/{stress_type}'
    names = ['overall', 'perSteps', 'own',]
    _gif_list(path_base, names, duration, loop)
        
def _img_to_gif_heights_1D(path, duration=500, loop=2):
    path_base = f'{path}/heights'
    names = ['heights_curve', 'ordered_heights_curve']
    _gif_list(path_base, names, duration, loop)
def _img_to_gif_heights_2D(path, duration=500, loop=2):
    path_base = f'{path}/heights'
    names = ['heights2D']
    _gif_list(path_base, names, duration, loop)
    
def _img_to_gif_energies(path, duration=500, loop=2):
    path_base = f'{path}/energies'
    names = ['energies']
    _gif_list(path_base, names, duration, loop)
    
def _write_rows(path, data):
    with open(path,"w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    
    
'''     
def img_to_gif(loop=2, fps=2):
    #fps: frames per seconds
    #loop: -1: no loop | 0: infinit loop | 1: loop once (see image twice) | 2: loop twice | ...

    # Stresses 2D
    path_base = f'{path}/{deployment}_deployment/stresses/{stress_type}'
    for name in ['all', 'perSteps', 'own']:
        img_name_i = f'"{path_base}/jpg/*{name}*.jpg"'
        gif_name_i = f'"{path_base}/jpg/gif/{name}.gif"'
        !ffmpeg  -loglevel panic -f image2 -r $fps -pattern_type glob -i $img_name_i -loop $loop $gif_name_i
        img_name_i = f'{path_base}/png/*{name}*.png'
        gif_name_i = f'{path_base}/png/gif/{name}.gif'
        !gifski --quiet -o $gif_name_i --fps $fps --repeat $loop --quality 100 $img_name_i
    # height 2D
    path_base = f'{path}/{deployment}_deployment/heights'
    img_name_i = f'"{path_base}/jpg/heights2D*.jpg"'
    gif_name_i = f'"{path_base}/jpg/gif/heights2D.gif"'
    !ffmpeg  -loglevel panic -f image2 -r $fps -pattern_type glob -i $img_name_i -loop $loop $gif_name_i
    img_name_i = f'{path_base}/png/heights2D*.png'
    gif_name_i = f'{path_base}/png/gif/heights2D.gif'
    !gifski --quiet -o $gif_name_i --fps $fps --repeat $loop --quality 100 $img_name_i

    # Stresses curves
    path_ = f'{path}/{deployment}_deployment/stresses/{stress_type}'
    for graph in ['stress_curve','ordered_stress_scatter', 'stress_scatter']:
        img_name_i = f'"{path_}/jpg/{graph}*.jpg"'
        gif_name_i = f'"{path_}/jpg/gif/{graph}.gif"'
        !ffmpeg  -loglevel panic -f image2 -r $fps -pattern_type glob -i $img_name_i -loop $loop $gif_name_i
        img_name_i = f'{path_}/png/{graph}*.png'
        gif_name_i = f'{path_}/png/gif/{graph}.gif'
        !gifski --quiet -o $gif_name_i --fps $fps --repeat $loop --quality 100 $img_name_i
    # Energies and heights curve
    path_ = f'{path}/{deployment}_deployment'
    for path0, name in  zip(['energies', 'heights'], ['energies', 'heights_curve']):
        path_g = path_ + f'/{path0}'
        img_name_i = f'"{path_g}/jpg/{name}*.jpg"'
        gif_name_i = f'"{path_g}/jpg/gif/{name}.gif"'
        !ffmpeg  -loglevel panic -f image2 -r $fps -pattern_type glob -i $img_name_i -loop $loop $gif_name_i
        img_name_i = f'{path_g}/png/{name}*.png'
        gif_name_i = f'{path_g}/png/gif/{name}.gif'
        !gifski --quiet -o $gif_name_i --fps $fps --repeat $loop --quality 100 $img_name_i
'''