from helpers_images import read_metadata
from deployment import regular_grid

def deploy_inputs_grid(path, deployments):
    datas = _read_inputs(path)
    
    category = 'regular_grid'
    for name, degree, rows, cols, steps, active_cells, target_percents in datas:
        print(f'name: {name}, ({rows=}, {cols=}):')
        for deployment in deployments:
            regular_grid(degree, rows, cols, category, name, steps, deployment, active_cells, target_percents)
            print(f'      done with {deployment=}')
            
def _read_inputs(path):
    lines_per_struct = 8
    with open(path, 'r') as f:
        nb_lines = len(f.read().splitlines())
        
    datas = []
    for i in range(0, nb_lines, lines_per_struct):
        datas.append(list(read_metadata(path, start_line=i)))
    return datas