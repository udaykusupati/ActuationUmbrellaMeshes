from tools import *
from images import *
from figure_2D import plot_undeployed_2D, projection2D
from RegularGrid import RegularGrid

import mesh
import sys
sys.path.append('..')
import configuration
sys.path.append('../UmbrellaGen')
import grid_gen

giff_total_time = 5000

def regular_grid(degree, rows, cols, category, name, steps, deployment, active_cells, target_percents,
                 heights_fct= None, min_height=64, verbose=False):
    folder_name , path = create_dir_hierarchy(category,
                                              degree,
                                              rows,
                                              cols,
                                              deployment,
                                              name)
    
    grid = RegularGrid(degree=degree,
                       rows=rows,
                       cols=cols,
                       height_fct=heights_fct,
                       min_height=min_height)
    
    grid.generate_mesh(folder_name, verbose=verbose)
    
    _deploy(path, folder_name, grid.input_data, grid.curr_um, degree, rows, cols, steps,
            active_cells, target_percents, grid.init_heights, grid.plate_thickness, deployment,
            verbose)

def non_regular_grid(mesh_path, degree, category, name, steps, deployment, active_cells, target_percents,
                     heights_fct=None, min_height=64, verbose=False):
    
    # default value for non-regular grid
    rows=cols=0
    
    folder_name , path = create_dir_hierarchy(category,
                                              degree,
                                              rows,
                                              cols,
                                              deployment,
                                              name)
    
    if mesh_path.endswith('.obj'):
        # read existing mesh
        base_mesh = mesh.Mesh(mesh_path)
        V, F = base_mesh.vertices(), base_mesh.elements()
        V_3d = np.zeros((len(V), 3))
        V_3d[:, :2] = V
        edge_length = np.linalg.norm(V_3d[F][0, 0] - V_3d[F][0, 1])
        numUmbrellas = len(F)
        base_mesh_gen = V_3d.tolist(), F.tolist()

        grid_gen.genUmbrellaWithHeights(degree, rows, cols,
                                        height_scales=None if heights_fct==None else heights_fct(numUmbrellas),
                                        minHeight=min_height,
                                        base_mesh=base_mesh_gen,
                                        edge_length=edge_length,
                                        json_filename=folder_name)
    
        input_path = '../UmbrellaGen/{}.json.gz'.format(folder_name)
    elif mesh_path.endswith('.json.gz'): input_path = mesh_path
    else: raise ValueError(f'mesh file not recognized: {mesh_path}, the extention should either be `.obj` or `.json.gz`.')
    
    io, input_data, target_mesh, curr_um, plate_thickness_scaled, target_height_multiplier = \
        configuration.parse_input(input_path, handleBoundary = False, isHex = (degree == 6), use_target_surface = False)

    init_heights = curr_um.umbrellaHeights
    
    _deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
            active_cells, target_percents, init_heights, plate_thickness_scaled, deployment,
            verbose)


def _deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
            active_cells, target_percents, init_heights, plate_thickness, deployment,
            verbose=False):
    
    write_metadata(path, folder_name, degree, rows, cols, steps, active_cells, target_percents)
    
    for phase, (active_c, target_p) in enumerate(zip(active_cells, target_percents)):
        plot_undeployed_2D(input_data,
                           curr_um,
                           show_height=False,
                           active_cells=active_c,
                           target_percents=target_p,
                           file_name = path+f'/phase{phase+1}_undeployed.png',
                           show_plot=verbose)
        '''
        # should be done after deployement, but required for each phase...
        # -not so meaningful-
        if (rows==1 or cols==1):
            projection2D(input_data['umbrella_connectivity'],
                         curr_um, active_c, target_p,
                         file_name=path+f'/phase{phase+1}_projection2D.png', show_plot = False)
        '''
    if verbose: print(f'Deployment')
    deploy_in_phase(curr_um,
                    input_data['umbrella_connectivity'],
                    init_heights,
                    plate_thickness,
                    active_cells,
                    target_percents,
                    deployment,
                    path,
                    steps=steps,
                    verbose=verbose)
    
    img_duration = giff_total_time/steps
    stress_types = ['VonMises','maxBending','Twisting']

    nb_phases = len(active_cells)
    for stress_type in stress_types:
        if verbose: print(f'\n-> generate images for {stress_type}.')
        for phase in range(nb_phases):
            if verbose: print(f'  phase: {phase:0>2}')
            generate_stresses_2D(path,
                                 phase,
                                 deployment,
                                 stress_type=stress_type,
                                 show_percent=False,
                                 show_plot=False,
                                 verbose=verbose)
        if verbose: print(f'  generate 1D images:')
        generate_stresses_1D([path],
                             [deployment],
                             stress_type=stress_type,
                             verbose=verbose)
        
    if verbose: print(f'  generate images for heights and energies.')
    for phase in range(nb_phases):
        if verbose: print(f'phase: {phase:0>2}')
        generate_heights_2D(path, phase, deployment,
                            verbose=verbose)
    if verbose: print(f'  generate 1D images:')
    generate_heightsEnergies_1D([path],
                                [deployment],
                                verbose=verbose)
    
    if verbose: print(f'\n-> generate GIFs.')
    img_to_gif(path, deployment, stress_type, duration=img_duration, loop=2, verbose=verbose)

