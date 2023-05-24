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

# ======================================================================
# =============================================================== GRID =
# ======================================================================
def set_actives_dep_weights(numUmbrellas, *active_cells,
                            dep_factors = np.logspace(-4, 0, 5)):
    weights_per_cell = np.zeros(numUmbrellas)
    weights_per_cell[active_cells] = 1
    return np.einsum('i, j -> ij', dep_factors, weights_per_cell)

def set_target_height(numUmbrellas, thickness, active_cell, target_heights):
    target_height_multiplier = np.ones(numUmbrellas)*-1
    target_height_multiplier[active_cell] = [h/thickness for h in target_heights]
    return target_height_multiplier

def percent_to_height(init_height, thickness, indexes, percents):
    return [(1-percent/100)*(init_height[idx]-thickness)+thickness
            for percent,idx in zip(percents,indexes)]

def height_to_percent(init_height, thickness, indexes, heights):
    return [max(0, min((1-(height-thickness)/(init_height[idx]-thickness))*100, 100))
            for height,idx in zip(heights[indexes], indexes)]

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

    os.makedirs(path+f'/undeployed')
    
    for type_ in ['heights', 'energies']:
        _sub_folder(path_dep+f'/{type_}')

    for s_type in help_.get_stress_types():
        _sub_folder(path_dep+f'/stresses/{s_type}')
    
    return folder_name, path

def write_metadata(path, folder_name, degree, rows, cols, steps, active_cells, target_percents):
    with open(path+"/metadata.txt", "w") as f:
        f.write("Name  : " + folder_name + '\n')
        f.write("Degree: " + str(degree) + '\n')
        f.write("Rows  : " + str(rows)   + '\n')
        f.write("Cols  : " + str(cols)   + '\n')
        f.write("Steps : " + str(steps)      + '\n')
        f.write("Phases: " + str(len(active_cells)) + '\n')
        f.write("Active Cells   : " + str(active_cells)    + '\n')
        f.write("Target Percents: " + str(target_percents) + '\n')
    

def deploy_in_phase(curr_um, connectivity, init_heights, plate_thickness, active_cells_all, target_percents_all, deployment, path,
                    steps=[10], verbose=True):

    # write Connectivity and Undeployed Units positions
    _write_rows(path+'/connectivity.csv', connectivity)
    _write_rows(path+'/position.csv', get_center_position(curr_um))

    prev_percents = [0]*len(target_percents_all[0]) # start with undeployed units
    max_phase = len(active_cells_all)-1

    stress_types = help_.get_stress_types()

    # to record max/min stresses during deployment:
    max_stresses_all = {}
    for stress_type in stress_types:
        max_stresses_all[stress_type] = np.zeros((curr_um.numUmbrellas(), curr_um.numUmbrellas()))

    nb_phases = len(active_cells_all)
    for phase, (active_cells, target_percents) in enumerate(zip(active_cells_all, target_percents_all)):
        if verbose:
            print(f'\n==== PHASE {phase+1:0>2}/{nb_phases:0>2} ====')
            print(f'{active_cells    = }')
            print(f'{target_percents = }')
        
        dep_weights = set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
        
        # to record values to write
        percents_per_steps = []
        heights            = []
        energies           = []

        # header for energies.csv
        energies.append(allEnergies(curr_um).keys()) # file will be read back as dict

        for s in range(steps[phase]+1):
            target_percents_s = _get_target_percents_s(deployment, s, steps[phase], prev_percents, target_percents)
            percents_per_steps.append(target_percents_s)

            target_heights = percent_to_height(init_heights, plate_thickness, active_cells, target_percents_s)
            target_height_multiplier = set_target_height(curr_um.numUmbrellas(), plate_thickness, active_cells, target_heights)
            success, _ = deploy_umbrella_pin_rigid_motion(curr_um,
                                                          plate_thickness,
                                                          target_height_multiplier,
                                                          dep_weights=dep_weights)
            if success:
                heights.append(curr_um.umbrellaHeights)
                energies.append(list(allEnergies(curr_um).values()))
                for stress_type in stress_types:
                    path_stresses = path+f'/{deployment}_deployment/stresses/{stress_type}/values'
                    max_stress_step = help_.get_max_stress_matrix(curr_um, stress_type)
                    # write max Stresses for this step
                    _write_rows(path_stresses+f'/phase{phase+1:0>2}_step{s:0>2}.csv',
                                max_stress_step)
                    
                    # retain max deployment stresses
                    is_max = max_stress_step > max_stresses_all[stress_type]
                    max_stresses_all[stress_type][is_max] = max_stress_step[is_max]

                if verbose: print(f'step {s:0>2}/{steps[phase]:0>2} saved.')

            else: raise ValueError(f'did not converge at step {s}.')

        # write Heights and Energies
        _write_rows(path+f'/{deployment}_deployment/heights/values/phase{phase+1:0>2}_heights.csv',
                    np.array(heights)-plate_thickness)
        _write_rows(path+f'/{deployment}_deployment/energies/values/phase{phase+1:0>2}_energies.csv',
                    energies)
        # write Percnets and Active Units
        with open(path+f'/{deployment}_deployment/heights/values/phase{phase+1:0>2}_percents.csv',"w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(active_cells)
            writer.writerows(percents_per_steps)
        
        # get openning percent of next phase's active units
        if phase < max_phase:
            prev_percents = height_to_percent(init_heights, plate_thickness, active_cells_all[phase+1], curr_um.umbrellaHeights)
    
    # write max deployment stress
    for stress_type in stress_types:
        path_stresses = path+f'/{deployment}_deployment/stresses/{stress_type}/values'
        _write_rows(path_stresses+f'/max_stresses.csv', max_stresses_all[stress_type])


        
def gif_to_img_duration(gif_duration, steps):
    return 1000*gif_duration/(steps+1)
        
def img_to_gif(path, deployment, stress_type, steps, duration=500, loop=0, verbose=False):
    img_to_gif_undeployed(f'{path}/undeployed', steps, duration=duration, loop=loop)
    
    path = f'{path}/{deployment}_deployment'
    img_to_gif_1D(path, stress_type, duration=duration, loop=loop)
    if verbose: print('gif 1D done')
    img_to_gif_2D(path, stress_type, duration=duration, loop=loop)
    if verbose: print('gif 2D done')

def img_to_gif_1D(path, stress_type, duration=500, loop=0):
    img_to_gif_stress_1D (path, stress_type, duration=duration, loop=loop)
    img_to_gif_heights_1D(path,              duration=duration, loop=loop)
    img_to_gif_energies  (path,              duration=duration, loop=loop)
    
def img_to_gif_2D(path, stress_type, duration=500, loop=2):
    _img_to_gif_stress_2D (path, stress_type, duration=duration, loop=loop)
    _img_to_gif_heights_2D(path,              duration=duration, loop=loop)
    
def img_to_gif_all(path, gif_name, duration=500, loop=2, steps=None):
    '''
    create gif `gif_name` with all files at path
    '''
    _create_gif(f'{path}/*.*', f'{path}/{gif_name}.gif', duration, loop, steps)

def img_to_gif_stress_1D(path, stress_type, duration=500, loop=0):
    path_base = f'{path}/stresses/{stress_type}'
    names = ['stress_curve','ordered_stress_scatter', 'stress_scatter']
    _gif_list(path_base, names, duration, loop)
    
def img_to_gif_heights_1D(path, duration=500, loop=0):
    path_base = f'{path}/heights'
    names = ['heights_curve', 'ordered_heights_curve']
    _gif_list(path_base, names, duration, loop)

def img_to_gif_energies(path, duration=500, loop=0):
    path_base = f'{path}/energies'
    names = ['energies']
    _gif_list(path_base, names, duration, loop)
    
def img_to_gif_undeployed(path, steps, duration=500, loop=0):
    img_to_gif_all(path, 'all_undeployed', duration, loop, steps)

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

# ---------------------------------------------------------------------- helpers

def _sub_folder(path):
    os.makedirs(path+'/values')
    os.makedirs(path+'/png/gif')
    os.makedirs(path+'/jpg/gif')

def _get_target_percents_s(deployment, s, steps, prev_target_percents, target_percents):
    if   deployment=='linear':
        # use of max/min because of numerical/deployment uncertainties
        return [max(0, min(prev_p+(p-prev_p)*s/steps, 100)) for prev_p, p in zip(prev_target_percents,target_percents)]
    else : raise ValueError(f'deployment unknown: {deployment} (should be:\'linear\')')

    # following deployment strategies need to be adapted for multi-phase deployment
    '''
    # no more interesting from deploymeny in phases
    elif deployment=='incremental':
        return [min(p, 100*s/steps) for p in target_percents]
    elif deployment=='max':
        max_p = max(target_percents)
        return [min(p, max_p*s/steps) for p in target_percents]
    else : raise ValueError(f'deployment unknown: {deployment} (choises:\'linear\',\'incremental\',\'max\')')
    '''

def _create_gif(images_name, gif_name, duration, loop, steps=None):
    factor = 1 if steps==None else (sum(steps)+1)/len(steps) # for undeployed meshes
    frames = [Image.open(image) for image in sorted(glob.glob(images_name))]
    frame_one = frames[0]
    frame_one.save(gif_name, append_images=frames[1:], save_all=True, duration=factor*duration, loop=loop)
    # duration: display time of each frame in milliseconds
    # loop : number of loops (0->infint, -1->no loop)
    
    
def _gif_list(path, ls, duration, loop):
    for name in ls:
        _create_gif(f'{path}/jpg/phase??_{name}*.jpg',
                   f'{path}/jpg/gif/{name}.gif',
                   duration, loop)

def _img_to_gif_stress_2D(path, stress_type, duration=500, loop=0):
    path_base = f'{path}/stresses/{stress_type}'
    names = ['overall', 'perSteps', 'own',]
    _gif_list(path_base, names, duration, loop)
    
def _img_to_gif_heights_2D(path, duration=500, loop=0):
    path_base = f'{path}/heights'
    names = ['heights2D']
    _gif_list(path_base, names, duration, loop)

def _write_rows(path, data):
    with open(path,"w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    
    
''' 
should have gifski and ffmpeg installed

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