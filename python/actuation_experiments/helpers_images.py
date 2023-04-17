
import numpy as np
import csv

def read_metadata(path):
    with open(path+"/metadata.txt") as f:
        metadata = f.read().splitlines()
    degree = int(metadata[0][8:])
    rows = int(metadata[1][8:])
    cols = int(metadata[2][8:])
    deployment = metadata[3][13:].strip()
    steps = int(metadata[4][13:])
    active_cells = eval(metadata[5][18:])
    target_percents = eval(metadata[6][18:])
    return degree, rows, cols, deployment, steps, active_cells, target_percents

def read_results(path, deployment, stress_type):
    path_dep = path+f'/{deployment}_deployment'
    path_stresses = path_dep+f'/{stress_type}/values'

    connectivity  = [[int(i), int(j)] for [i,j] in _read_csv(path+'/connectivity.csv')]
    init_position = _read_csv(path+'/position.csv')

    heights  = _read_csv(path_dep+'/heights.csv')
    actuation = _read_csv(path_dep+'/percents.csv')
    # format data:
    active_cells = [int(c) for c in actuation[0]] # 1st line is the activated cells indexes[0]
    percents_per_steps = actuation[1:]

    stresses = []
    for i in range(len(heights)):
        stresses.append(_read_csv(path_stresses+f'/step_{i:0>2}.csv'))

    return connectivity,\
           np.array(init_position),\
           heights,\
           active_cells,\
           percents_per_steps,\
           np.array(stresses)

def get_indexes(degree, rows, cols):
    if degree==3:
        indexes = []
        for r in range(rows):
            for c in range(r*2*cols, (r+1)*2*cols, 2):
                if r%2 == 0:
                    indexes.extend([c, c+1])
                else:
                    indexes.extend([c+1, c])
        return indexes
    
    if degree==4:
        return list(range(rows*cols))
    

# ----------------------------------------------------------------------
# ------------------------------------------------------------ helpers -
# ----------------------------------------------------------------------
def _read_csv(path):
    out = []
    with open(path,"r", newline='') as csvfile:
        reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            out.append(row)
    return out