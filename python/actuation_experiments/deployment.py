from tools import *
from images import *
from figure_2D import plot2D, projection2D
from RegularGrid import RegularGrid

import mesh
import sys
sys.path.append('..')
import configuration
sys.path.append('../UmbrellaGen')
import grid_gen


def regular_grid(degree, rows, cols, category, name, steps, deployment, active_cells, target_percents,
                 heights_fct= None, min_height=64, verbose=False):
    folder_name , path = create_dir_hierarchy(category,
                                              degree,
                                              rows,
                                              cols,
                                              deployment,
                                              name)

    write_metadata(path, degree, rows, cols, steps, active_cells, target_percents)
    
    
    grid = RegularGrid(degree=degree,
                       rows=rows,
                       cols=cols,
                       height_fct=heights_fct,
                       min_height=min_height)
    
    grid.generate_mesh(folder_name, verbose=verbose)
    
    plot2D(grid.input_data,
           grid.curr_um,
           show_height=False,
           active_cells=active_cells,
           target_percents=target_percents,
           file_name = path+'/undeployed.png',
           show_plot=verbose)
    
    deploy_in_steps(grid.curr_um,
                grid.input_data,
                grid.init_heights,
                grid.plate_thickness,
                active_cells,
                target_percents,
                deployment,
                path,
                steps=steps,
                verbose=verbose)
    
    if (grid.rows==1 or grid.cols==1):
        projection2D(grid.input_data['umbrella_connectivity'],
                        grid.curr_um, active_cells, target_percents,
                        file_name=path+'/projection2D.png', show_plot = False)
    
    stress_type = 'maxBending'
    generate_2D(path,
                deployment,
                stress_type=stress_type,
                show_percent=False,
                show_plot=False,
                verbose=verbose)
    
    generate_1D([path],
                [deployment],
                save_dir = '', 
                stress_type=stress_type,
                show_percent=False,
                show_plot=False,
                verbose=verbose)
    
    img_to_gif(path, deployment, stress_type, duration=500, loop=2, verbose=verbose)
    
    return

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

    init_center_pos = get_center_position(curr_um)
    init_heights = curr_um.umbrellaHeights
    
    write_metadata(path, degree, rows, cols, steps, active_cells, target_percents)
    
    plot2D(input_data,
           curr_um,
           show_height=False,
           active_cells=active_cells,
           target_percents=target_percents,
           file_name = path+'/undeployed.png',
           show_plot=verbose)
    
    dep_weights              = set_actives_dep_weights(curr_um.numUmbrellas(), active_cells)
    target_heights           = percent_to_height(init_heights, plate_thickness_scaled, active_cells, target_percents)
    target_height_multiplier = set_target_height(curr_um.numUmbrellas(), active_cells, target_heights)
    
    deploy_in_steps(curr_um,
                    input_data,
                    init_heights,
                    plate_thickness_scaled,
                    active_cells,
                    target_percents,
                    deployment,
                    path,
                    steps=steps,
                    verbose=verbose)
    
    stress_type='maxBending'
    generate_2D(path,
                deployment,
                stress_type=stress_type,
                show_percent=False,
                show_plot=False,
                verbose=verbose)
    generate_1D([path], [deployment],
                save_dir = '',
                stress_type=stress_type,
                show_percent=False,
                show_plot=False,
                verbose=verbose)
    img_to_gif(path, deployment, stress_type, duration=500, loop=2, verbose=verbose)
    
    
def read_inputs(path):
    with open(path) as f:
        inputs = f.read().splitlines()
    name = inputs[0][6:]
    degree = int(inputs[1][8:])
    rows = int(inputs[2][6:])
    cols = int(inputs[3][6:])
    steps = int(inputs[4][7:])
    active_cells = eval(inputs[5][14:])
    target_percents = eval(inputs[6][17:])
    return name, degree, rows, cols, steps, active_cells, target_percents