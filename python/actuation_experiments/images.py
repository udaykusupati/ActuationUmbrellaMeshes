import os
import numpy as np
import matplotlib.pyplot as plt

import helpers_images
import figure_2D 

# ======================================================================
# ==================================================== GENERATE IMAGES =
# ======================================================================
def generate_2D(path, deployment,
                stress_type='maxBending', show_percent=False, show_plot=False, verbose=False):
    
    connectivity,        \
    init_center_pos,     \
    heights,             \
    active_cells,        \
    percents_per_steps,  \
    stresses_per_steps,  \
    el_energies =        \
        helpers_images.read_results(path, deployment, stress_type)

    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()

    c = np.array(connectivity)
    max_stress_per_arm = stresses_per_steps.transpose()[c[:,0], c[:,1]].max(axis=1)
    
    title_s = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
    title_h = '{{}} - {:0>3.0f}% deployed'
    path_stresses = []
    path_heights = []
    for f in ['jpg', 'png']:
        path_stresses.append(f'{path}/{deployment}_deployment/stresses/{stress_type}/{f}'+'/{{}}_{:0>3.0f}Deployed'+f'.{f}')
        path_heights.append(f'{path}/{deployment}_deployment/heights/{f}'+'/{{}}_{:0>3.0f}Deployed'+f'.{f}')
    
    # random perturbations does affect undeployed state
    deployed = False
    
    for step, (s_matrix, percents) in enumerate(zip(stresses_per_steps,percents_per_steps)):
        min_stress_step = s_matrix[s_matrix != 0].min()
        max_stress_step = s_matrix[s_matrix != 0].max()
        max_stresses.append(max_stress_step)

        title_step_s = title_s.format(step/steps*100, max_stress_step)
        title_step_h = title_h.format(step/steps*100, max_stress_step)
        path_s = []
        path_h = []
        for s, h in zip(path_stresses, path_heights):
            path_s.append(s.format(step/steps*100))
            path_h.append(h.format(step/steps*100))
        
        # normalized with general extrems values
        figure_2D.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_all, deployed*max_stress_all, show_plot, title_step_s.format('all'), path_s, 'all')
        
        # normalized with step extrems values
        figure_2D.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_step, deployed*max_stress_step, show_plot, title_step_s.format('step'), path_s, 'perSteps')

        # normalized with own extrems values
        figure_2D.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, 0, deployed*max_stress_per_arm, show_plot, title_step_s.format('own'), path_s, 'own')
        
        # heights
        figure_2D.fig_height2D(connectivity, active_cells, heights[step], init_center_pos, show_plot, title_step_h.format('heights'), path_h, 'heights2D' )
        
        # perturbations do not affect deployed state
        deployed = True

        if verbose: print(f'2D images for step {step:0>2} successfully saved.')

def generate_1D(paths, deployments,
                save_dir = '', stress_type='maxBending', show_percent=False,
                show_plot=False, verbose=False):
    

    degrees      = []
    rows         = []
    cols         = []
    heights      = []
    indexes      = []
    stresses     = []
    energies     = []
    active_cells = []
    percents_per_steps = []

    for path, dep in zip(paths, deployments):
            degree_, rows_, cols_, *_  = helpers_images.read_metadata(path)
            connectivity_,_, heights_, active_cells_, percents_per_steps_, stresses_,  el_energies_ =\
                    helpers_images.read_results(path, dep, stress_type)

            degrees.append(degree_)
            rows.append(rows_)
            cols.append(cols_)
            heights.append(heights_)
            indexes.append(helpers_images.get_indexes(degree_, rows_, cols_))
            c = np.array(connectivity_)
            stresses.append([s[c[:,0], c[:,1]]for s in stresses_])
            energies.append(el_energies_)
            active_cells.append(active_cells_)
            percents_per_steps.append(percents_per_steps_)

    if len(paths) == 1:
         save_dir = paths[0] + f'/{deployments[0]}_deployment'
    if save_dir!= '':
        if not os.path.exists(save_dir):
              os.makedirs(save_dir+ '/energies/jpg/gif')
              os.makedirs(save_dir+ '/energies/png/gif')
              os.makedirs(save_dir+ '/heights/jpg/gif')
              os.makedirs(save_dir+ '/heights/png/gif')
              os.makedirs(save_dir+ f'/stresses/{stress_type}/jpg/gif')
              os.makedirs(save_dir+ f'/stresses/{stress_type}/png/gif')
        path_e = save_dir + '/energies'
        path_h = save_dir + '/heights'
        path_s = save_dir + f'/stresses/{stress_type}'
    else : path_e = path_h = path_s = ''

    figure_2D.figs_heights_curve(indexes,
                                 heights,
                                 active_cells,
                                 percents_per_steps,
                                 _add_jpg_png(path_h),
                                 show_active = True,
                                 ordered=False,
                                 show_plot = show_plot)
    if verbose: print('heights curve done')
    
    figure_2D.figs_heights_curve(indexes,
                                 heights,
                                 active_cells,
                                 percents_per_steps,
                                 _add_jpg_png(path_h),
                                 show_active = True,
                                 ordered=True,
                                 show_plot = show_plot)
    if verbose: print('ordered heights curve done')
    
    figure_2D.figs_stress_curve(stresses, _add_jpg_png(path_s), show_plot = show_plot)
    if verbose: print('stress curve done')
    
    figure_2D.figs_stress_scatter(stresses, _add_jpg_png(path_s), show_plot = show_plot)
    if verbose: print('stress scatter done')
    
    figure_2D.figs_stress_scatter(stresses, _add_jpg_png(path_s), ordered=False, show_plot = show_plot)
    if verbose: print('ordered stress scatter done')
    
    figure_2D.figs_energy_curve(energies, _add_jpg_png(path_e), show_plot = show_plot)
    if verbose: print('energy curve done')


def _add_jpg_png(path=''):
    paths = []
    if path != '':
        for f in ['jpg', 'png']:
            paths.append(f'{path}/{f}'+'/{{}}_{:0>3.0f}Deployed'+f'.{f}')
    return paths
