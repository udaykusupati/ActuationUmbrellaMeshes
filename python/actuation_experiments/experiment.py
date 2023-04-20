from tools import *
from images import *
from figure_2D import plot2D
from RegularGrid import RegularGrid


def run(degree, rows, cols, category, name, steps, deployment, active_cells, target_percents, verbose=False):
    folder_name , path = \
    create_dir_hierarchy(category,
                         degree,
                         rows,
                         cols,
                         deployment,
                         name)

    write_metadata(path, degree, rows, cols, steps, active_cells, target_percents)
    
    
    grid = RegularGrid(degree, rows, cols)
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
                path,
                deployment,
                steps=steps,
                verbose=verbose,
                dep=deployment,
                prj2D=(grid.rows==1 or grid.cols==1))
    
    stress_type = 'maxBending'
    generate_2D(path,
            deployment,
            stress_type=stress_type,
            show_percent=False,
            show_plot=False,
            verbose=verbose)
    
    generate_1D([path],
            [deployment],
            stress_type=stress_type,
            show_percent=False,
            show_plot=False,
            verbose=verbose)
    
    img_to_gif(path, deployment, stress_type, duration=500, loop=2, verbose=verbose)
    
    return

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