import numpy as np
import matplotlib.pyplot as plt

import helpers_images
import helpers_plots

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
    stresses_per_steps = \
        helpers_images.read_results(path, deployment, stress_type)

    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()

    c = np.array(connectivity)
    max_stress_per_arm = stresses_per_steps.transpose()[c[:,0], c[:,1]].max(axis=1)
    
    title = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
    path_names = []
    for f in ['jpg', 'png']:
        path_names.append(f'{path}/{deployment}_deployment/{stress_type}/{f}'+'/{{}}_{:0>3.0f}Deployed'+f'.{f}')
    
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
        helpers_plots.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_all, deployed*max_stress_all, show_plot, title_s.format('all'), path_names_s, 'all')
        
        # normalized with step extrems values
        helpers_plots.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, deployed*min_stress_step, deployed*max_stress_step, show_plot, title_s.format('step'), path_names_s, 'perSteps')

        # normalized with own extrems values
        helpers_plots.fig_arm_stresses(connectivity, active_cells, percents, init_center_pos, show_percent,
                              s_matrix, 0, deployed*max_stress_per_arm, show_plot, title_s.format('own'), path_names_s, 'own')
        
        # perturbations do not affect deployed state
        deployed = True

        if verbose: print(f'Images for step {s:0>2} successfully saved.')

def generate_1D(paths, deployments,
                save_dir = '', stress_type='maxBending', show_percent=False, show_plot=False):
    

    degrees = []
    rows    = []
    cols    = []
    heights = []
    indexes = []
    stresses = []
    percents = []
    connectivity = []

    for path, dep in zip(paths, deployments):
        degree_, rows_, cols_, *_  = helpers_images.read_metadata(path)
        connectivity_,_, heights_, active_cells, percents_per_steps, stresses_ =\
              helpers_images.read_results(path, dep, stress_type)
        
        degrees.append(degree_)
        rows.append(rows_)
        cols.append(cols_)
        heights.append(heights_)
        indexes.append(helpers_images.get_indexes(degree_, rows_, cols_))
        c = np.array(connectivity_)
        stresses.append([s[c[:,0], c[:,1]]for s in stresses_])

    if save_dir == '' : save_dir = paths[0]

    helpers_plots.figs_heights(indexes, heights)
    helpers_plots.figs_stress_curve(max_stresses)
    helpers_plots.figs_stress_scatter(stresses)
    helpers_plots.figs_stress_scatter(stresses, ordered=False)
    # add red/green dot for heights
    # add plot for energies
    # reverse height for poercent deployment ?
    # save fgis





    steps = stresses_per_steps.shape[0]-1
    max_stresses = []
    stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
    min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()
    max_x, max_y = steps, 1.05*max_stress_all

    c = np.array(connectivity)
    max_stress_per_arm = stresses_per_steps.transpose()[c[:,0], c[:,1]].max(axis=1)
    
    title = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
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

        # ordered stresses
        help_.fig_stress_scatter(connectivity, deployed*s_matrix, max_y, show_plot, title_s.format('ordered stress'), path_names_s, 'ordered_scatter')
        
        # stress curve
        help_.fig_stress_curve(deployed*max_stresses, max_x, max_y, show_plot, title_s.format('deployment'), path_names_s, 'sPlot')
        
        # perturbations do not affect deployed state
        deployed = True





# def generate_1D(connectivity, active_cells, percents_per_steps, init_center_pos, stresses_per_steps, dir_name='./00_test',
#                  stress_type='maxBending', show_percent=False, show_plot=False):
#     steps = stresses_per_steps.shape[0]-1
#     max_stresses = []
#     stresses_all_nz = stresses_per_steps[stresses_per_steps != 0]
#     min_stress_all, max_stress_all = stresses_all_nz.min(), stresses_all_nz.max()
#     max_x, max_y = steps, 1.05*max_stress_all
    
#     title = '{{}} - {:0>3.0f}% deployed ('+f'{stress_type}:'+' {:0>6.2f})'
#     path_names = []
#     if dir_name == '': img_format = []
#     else:              img_format = ['jpg', 'png']
#     for f in img_format:
#         path_names.append(f'{dir_name}/{f}/{stress_type}_'+'{{}}_{:0>3.0f}Deployment'+f'.{f}')
    
#     # random perturbations does affect undeployed state
#     deployed = False
    
#     for s, (s_matrix, percents) in enumerate(zip(stresses_per_steps,percents_per_steps)):
#         max_stress_step = s_matrix[s_matrix != 0].max()

#         title_s = title.format(s/steps*100, max_stress_step)
#         path_names_s = []
#         for path in path_names:
#             path_names_s.append(path.format(s/steps*100))
        
#         help_.fig_stress_scatter(connectivity, deployed*s_matrix, max_y, show_plot, title_s.format('stress'), path_names_s, 'scatter', ordered=False)
        
#         # perturbations do not affect deployed state
#         deployed = True
