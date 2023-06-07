import os
import numpy as np
import matplotlib.pyplot as plt

import helpers_images
import figure_2D 

# ======================================================================
# ==================================================== GENERATE IMAGES =
# ======================================================================
def generate_stresses_2D(path, phase, deployment,
                stress_type='VonMises', show_percent=False,
                show_plot=False, verbose=False):
    # ensure plt.rcParams.update() to update for next plots
    figure_2D.fig_empty()
    
    connectivity,        \
    init_center_pos,     \
    _,                   \
    active_cells,        \
    percents_per_steps,  \
    stresses_per_steps,  \
    max_stresses,        \
    _                    \
        = helpers_images.read_results(path, phase+1, deployment, stress_type)

    steps = stresses_per_steps.shape[0]-1
    max_stress_per_arm = max_stresses.transpose()[connectivity[:,0], connectivity[:,1]]
    max_stress_all = max_stress_per_arm.max()
    
    title_s = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
    path_stresses = []
    for f in ['jpg', 'png']:
        path_stresses.append(f'{path}/{deployment}_deployment/stresses/{stress_type}/{f}/phase{(phase+1):0>2}_'+'{{}}_{:0>3.0f}Deployed'+f'.{f}')
    
    # ensure stress 0 for undeployed state (random perturbation)
    deployed = False if phase==0 else True # only undeployed at phase 0, step 0
    
    for step, (s_matrix, percents) in enumerate(zip(stresses_per_steps, percents_per_steps)):
        if phase>0 and step==0: continue
        min_stress_step = s_matrix.min()
        max_stress_step = s_matrix.max()
        # min_stress_step = s_matrix[s_matrix!=0].min()
        # max_stress_step = s_matrix[s_matrix!=0].max()
        
        # avoid non null stress for non deployed
        if phase==0 and step==0:
            min_stress_step = 0.0
            max_stress_step = 0.0

        title_step_s = title_s.format(step/steps*100, max_stress_step)
        path_s = []
        for s in path_stresses:
            path_s.append(s.format(step/steps*100))
        
        # normalized with general extrems values
        figure_2D.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, 0, deployed*max_stress_all, show_plot, title_step_s.format('overall'), path_s, 'overall')
        
        # normalized with step extrems values
        figure_2D.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_step, deployed*max_stress_step, show_plot, title_step_s.format('step'), path_s, 'perSteps')

        # normalized with own extrems values
        figure_2D.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, 0, deployed*max_stress_per_arm, show_plot, title_step_s.format('own'), path_s, 'own')
         
        # perturbations do not affect deployed state
        deployed = True

        if verbose: print(f'2D stresses images for step {step:0>2} successfully saved.')


def generate_heights_2D(path, phase, deployment, plate_thickness, init_heights,
                        stress_type='VonMises', show_plot=False, verbose=False):
    # ensure plt.rcParams.update() to update for next plots
    figure_2D.fig_empty()
    
    connectivity,        \
    init_center_pos,     \
    heights,             \
    active_cells,        \
    _,                   \
    _,                   \
    _,                   \
    _                    \
        = helpers_images.read_results(path, phase+1, deployment, stress_type)

    steps = heights.shape[0]-1
    
    title_h = '{{}} - {:0>3.0f}% deployed'
    path_heights = []
    for f in ['jpg', 'png']:
        path_heights.append(f'{path}/{deployment}_deployment/heights/{f}/phase{(phase+1):0>2}_'+'{{}}_{:0>3.0f}Deployed'+f'.{f}')
    
    
    for step, heights_, in enumerate(heights):
        if phase>0 and step==0: continue
 
        title_step_h = title_h.format(step/steps*100)
        path_h = []
        for h in path_heights:
            path_h.append(h.format(step/steps*100))
        
        # heights
        figure_2D.fig_height2D(connectivity, active_cells, plate_thickness, init_heights, heights_, init_center_pos, show_plot, title_step_h.format('heights'), path_h, 'heights2D' )
        if verbose: print(f'2D heights images for step {step:0>2} successfully saved.')

# --------------------------------------------------------------------------------
# ----------------------------------------------------------------------------- 1D
# --------------------------------------------------------------------------------
def generate_stresses_1D(paths, deployments,
                         save_dir = '', stress_type='VonMises', show_percent=False,
                         show_plot=False, verbose=False):
    # ensure plt.rcParams.update() to update for next plots
    figure_2D.fig_empty()
    
    stresses    = []
    list_steps  = []
    list_phases = []
    
    nb_steps = 0
    nb_arms  = 0
    s_max    = 0
    for path, dep in zip(paths, deployments):
        _,_,_,_, steps, phases, *_  = helpers_images.read_metadata(path+'/metadata.txt')
        list_steps.append(steps)
        list_phases.append(phases)
        stresses_phase = []
        indexes_phase  = []
        if nb_steps == 0: nb_steps = sum(steps)
        for phase in range(phases):
            connectivity_,_,_,_,_, stresses_,max_stresses_,_ =\
                helpers_images.read_results(path, phase+1, dep, stress_type)
            if nb_arms == 0: nb_arms = connectivity_.shape[0]
            stresses_phase.append([s[connectivity_[:,0], connectivity_[:,1]]for s in stresses_])
            if max_stresses_.max() > s_max: s_max = max_stresses_.max()
        stresses.append(stresses_phase)
    
    max_y = s_max*1.05 # add some margin for plot
    max_x_curve   = nb_steps
    max_x_scatter = nb_arms
    
    if len(paths) == 1:
         save_dir = paths[0] + f'/{deployments[0]}_deployment'
    if save_dir!= '':
        if not os.path.exists(save_dir+ f'/stresses/{stress_type}'):
              os.makedirs(save_dir+ f'/stresses/{stress_type}/jpg/gif')
              os.makedirs(save_dir+ f'/stresses/{stress_type}/png/gif')
        path_s = save_dir + f'/stresses/{stress_type}'
    else : path_s = ''
    
    figure_2D.figs_stress_curve(stresses, (max_x_curve ,max_y), _add_jpg_png(path_s), show_plot = show_plot)
    if verbose: print('stress curve done')
    
    figure_2D.figs_stress_scatter(stresses, (max_x_scatter ,max_y), _add_jpg_png(path_s), show_plot = show_plot)
    if verbose: print('stress scatter done')
    
    figure_2D.figs_stress_scatter(stresses, (max_x_scatter ,max_y), _add_jpg_png(path_s), ordered=False, show_plot = show_plot)
    if verbose: print('ordered stress scatter done')
    

def generate_heightsEnergies_1D(paths, deployments,
                                save_dir = '', stress_type='VonMises', show_plot=False, verbose=False):
    # ensure plt.rcParams.update() to update for next plots
    figure_2D.fig_empty()
    
    heights            = []
    indexes            = []
    energies           = []
    active_cells       = []
    percents_per_steps = []
    
    h_max = 0
    e_max = 0
    nb_units = 0
    nb_steps = 0
    for path, dep in zip(paths, deployments):
        _, degree_, rows_, cols_, steps, phases, *_  = helpers_images.read_metadata(path+'/metadata.txt')
        heights_phase            = []
        indexes_phase            = []
        energies_phase           = []
        active_cells_phase       = []
        percents_per_steps_phase = []
        if nb_steps == 0: nb_steps = sum(steps)
        for phase in range(phases):
            _,_, heights_, active_cells_, percents_per_steps_, stresses_, _, el_energies_ =\
                    helpers_images.read_results(path, phase+1, dep, stress_type)
            if nb_units == 0: nb_units = stresses_.shape[1]
            heights_phase.append(heights_)
            if rows_==0 or cols_==0: indexes_phase.append(list(range(len(heights_[0])))) # external mesh
            else: indexes_phase.append(helpers_images.get_indexes(degree_, rows_, cols_))
            energies_phase.append(el_energies_)
            active_cells_phase.append(active_cells_)
            percents_per_steps_phase.append(percents_per_steps_)
            
            # max values
            if heights_.max() > h_max: h_max = heights_.max()
            if el_energies_.max() > e_max: e_max = el_energies_.max()
        heights.append(heights_phase)
        indexes.append(indexes_phase)
        energies.append(energies_phase)
        active_cells.append(active_cells_phase)
        percents_per_steps.append(percents_per_steps_phase)
        
    h_xylim = (nb_units, 1.05*h_max)
    e_xylim = (nb_steps, 1.05*e_max)
                
    if len(paths) == 1:
         save_dir = paths[0] + f'/{deployments[0]}_deployment'
    if save_dir!= '':
        if not os.path.exists(save_dir+ '/energies'):
              os.makedirs(save_dir+ '/energies/jpg/gif')
              os.makedirs(save_dir+ '/energies/png/gif')
              os.makedirs(save_dir+ '/heights/jpg/gif')
              os.makedirs(save_dir+ '/heights/png/gif')
        path_e = save_dir + '/energies'
        path_h = save_dir + '/heights'
    else : path_e = path_h = ''
    
    figure_2D.figs_heights_curve(indexes,
                                 heights,
                                 h_xylim,
                                 active_cells,
                                 percents_per_steps,
                                 _add_jpg_png(path_h),
                                 show_active = True,
                                 ordered=False,
                                 show_plot = show_plot)
    if verbose: print('heights curve done')
    
    figure_2D.figs_heights_curve(indexes,
                                 heights,
                                 h_xylim,
                                 active_cells,
                                 percents_per_steps,
                                 _add_jpg_png(path_h),
                                 show_active = True,
                                 ordered=True,
                                 show_plot = show_plot)
    if verbose: print('ordered heights curve done')
    
    figure_2D.figs_energy_curve(energies, e_xylim, _add_jpg_png(path_e), show_plot = show_plot)
    if verbose: print('energy curve done')

def _add_jpg_png(path=''):
    paths = []
    if path != '':
        for f in ['jpg', 'png']:
            paths.append(f'{path}/{f}/'+'phase{:0>2}_{{}}_{:0>3.0f}Deployed'+f'.{f}')
    return paths
