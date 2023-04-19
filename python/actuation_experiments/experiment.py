from tools import *
from figure_2D import *
from images import *
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