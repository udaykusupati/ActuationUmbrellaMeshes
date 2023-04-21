
import numpy as np
import csv

def read_metadata(path, start_line=0):
    with open(path) as f:
        metadata = f.read().splitlines()
    name            =      metadata[start_line+0][8:]
    degree          = int (metadata[start_line+1][8:])
    rows            = int (metadata[start_line+2][8:])
    cols            = int (metadata[start_line+3][8:])
    steps           = int (metadata[start_line+4][8:])
    active_cells    = eval(metadata[start_line+5][18:])
    target_percents = eval(metadata[start_line+6][18:])
    return name, degree, rows, cols, steps, active_cells, target_percents

def read_results(path, deployment, stress_type):
    path_dep = path+f'/{deployment}_deployment'
    path_stresses = path_dep+f'/stresses/{stress_type}/values'

    connectivity  = [[int(i), int(j)] for [i,j] in _read_csv(path+'/connectivity.csv')]
    init_position = _read_csv(path+'/position.csv')

    heights  = _read_csv(path_dep+'/heights/values/heights.csv')
    actuation = _read_csv(path_dep+'/heights/values/percents.csv')
    # format data:
    active_cells = [int(c) for c in actuation[0]] # 1st line is the activated cells indexes[0]
    percents_per_steps = actuation[1:]

    # return only elastic energy
    el_energies = _read_csv_dic(path_dep+'/energies/values/energies.csv')


    stresses = []
    for i in range(len(heights)):
        stresses.append(_read_csv(path_stresses+f'/step_{i:0>2}.csv'))

    return connectivity,\
           np.array(init_position),\
           heights,\
           active_cells,\
           percents_per_steps,\
           np.array(stresses),\
           el_energies

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

def _read_csv_dic(path):
    # # return a list of dict
    # return csv.DictReader(open(path))
    # keys: ['Full','Elastic','Deployment','Repulsion','Attraction','AngleBoundPenalty']
    
    # return only elastic energy (index 0):
    datas = []
    with open(path, "r") as csvfile: 
        reader = csv.reader(csvfile)
        header = next(reader, None)
        for row in reader:
            datas.append(row)
    return [float(energies[1]) for energies in datas]
