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
from figure_2D import projection2D

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

def write_metadata(path, degree, rows, cols, steps, active_cells, target_percents):
    with open(path+"/metadata.txt", "w") as f:
        f.write("Degree: " + str(degree) + '\n')
        f.write("Rows  : " + str(rows)   + '\n')
        f.write("Cols  : " + str(cols)   + '\n')
        f.write("Steps      : " + str(steps)      + '\n')
        f.write("Active Cells    : " + str(active_cells)    + '\n')
        f.write("Target Percents : " + str(target_percents) + '\n')

def deploy_in_steps(curr_um, input_data, init_heights, plate_thickness, active_cells, target_percents, path, deployment,
                    steps=10, verbose=True, dep='linear', prj2D=False):
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
    
    if prj2D: projection2D(input_data['umbrella_connectivity'],
                            curr_um, active_cells, target_percents,
                            file_name=path+'/projection2D.png', show_plot = False)
    
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

def img_to_gif(path, deployment, stress_type, duration=500, loop=2, verbose=False):
    img_to_gif_stress  (path, deployment, stress_type, duration, loop)
    if verbose: print('gif for stresses done')
    img_to_gif_heights (path, deployment, duration, loop)
    if verbose: print('gif for heights done')
    img_to_gif_energies(path, deployment, duration, loop)
    if verbose: print('gif for energy done')
    
def img_to_gif_stress(path, deployment, stress_type, duration=500, loop=2):
    path_base = f'{path}/{deployment}_deployment/stresses/{stress_type}'
    for name in ['all', 'perSteps', 'own', 'stress_curve','ordered_stress_scatter', 'stress_scatter']:
        create_gif(f'{path_base}/jpg/{name}*.jpg',
                   f'{path_base}/jpg/gif/{name}.gif',
                   duration, loop)
        
def img_to_gif_heights(path, deployment, duration=500, loop=2):
    path_base = f'{path}/{deployment}_deployment/heights'
    for name in ['heights2D', 'heights_curve', 'ordered_heights_curve']:
        create_gif(f'{path_base}/jpg/{name}*.jpg',
                   f'{path_base}/jpg/gif/{name}.gif',
                   duration, loop)
    
def img_to_gif_energies(path, deployment, duration=500, loop=2):
    path_base = f'{path}/{deployment}_deployment/energies'
    create_gif(f'{path_base}/jpg/energies*.jpg',
               f'{path_base}/jpg/gif/energies.gif',
               duration, loop)

        
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


def create_gif(images_name, gif_name, duration, loop):
    frames = [Image.open(image) for image in sorted(glob.glob(images_name))]
    frame_one = frames[0]
    frame_one.save(gif_name, append_images=frames[1:], save_all=True, duration=duration, loop=loop)
    
    
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