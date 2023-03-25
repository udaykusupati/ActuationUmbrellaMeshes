## ===============
## UTILS
## ===============
def get_num_cells(rows, cols):
    return rows*(cols*2)

def get_border(rows, cols):
    step = 2*cols
    bot_left  = 0
    bot_right = step-1
    top_left  = step*(rows-1)
    top_right = step*rows
    return  list(range(bot_right)) + \
            list(range(bot_right, top_right-1, step)) + \
            list(range(step, top_left, step)) + \
            list(range(top_left,top_right))

def get_center(rows, cols):
    if rows%2==0:
        a = (rows)*(cols-1)
        b = (rows+1)*cols
        return [a-1, a, b-1, b]
    else:
        return [rows*cols-1, rows*cols]

## ===============
## HEIGHTS
## ===============
def height_cst(rows, cols):
    return [1]*get_num_cells(rows, cols)

def height_Ux(rows, cols):
    numUmbrellas = get_num_cells(rows, cols)
    heights = height_cst(rows, cols)
    for uid in range(numUmbrellas):
        heights[uid] += (0.01 * uid**2 + 0.01 * (numUmbrellas - 1 - uid)**2)
    heights = [h/min(heights) for h in heights]
    return heights

def height_Ux_nbCell(numUmbrellas):
    heights = [1]*numUmbrellas
    for uid in range(numUmbrellas):
        heights[uid] += (0.01 * uid**2 + 0.01 * (numUmbrellas - 1 - uid)**2)
    heights = [h/min(heights) for h in heights]
    return heights


## ===============
## PLOT 2D - close Umbrella
## ===============
import numpy as np
import matplotlib.pyplot as plt

def plot2D(input_data, curr_um, show_height=False, active_cells=[], target_percents=[]):
    # get cells' center coordinates
    nb_cell = curr_um.numUmbrellas()
    center_position = np.zeros([nb_cell, 3])
    for i in range(nb_cell):
        top_idx = curr_um.getUmbrellaCenterJi(i, 0)
        center_position[i] = curr_um.joint(top_idx).position

    # get cells' top plate edges coordinate
    vertices = input_data['base_mesh_v']
    edge = np.roll(np.insert(input_data['base_mesh_f'], 0, input_data['base_mesh_f'][:,0], axis=1), -1, axis=1)

    # == actual plotting ==
    fig = plt.figure(figsize=(15, 15))
    # plot cells' edges
    for e in vertices[edge]:
        plt.plot(e[:, 0], e[:, 1], color="lightblue")
    # plot cell index at its center
    for i, [x,y,z] in enumerate(center_position):
        plt.annotate(f'{i}', (x,y), ha='center', color='black')
    if show_height:
        h_max, h_min = max(curr_um.umbrellaHeights), min(curr_um.umbrellaHeights)
        for [x,y,z], h in zip(center_position, curr_um.umbrellaHeights):
            if h_max != h_min :
                r = (h-h_min)/(h_max-h_min)
                g = 0.1
                b = 1-r
            else:
                r,g,b = 0,0,0
            plt.annotate(f'[{h:.2f}]', (x,y), textcoords='offset points', xytext=(0,-12), ha='center', c=(r,g,b))
    # overwrite index if actve_cell
    if active_cells != []:
        for i, p in zip(active_cells, target_percents):
            g = p/100
            r = 1-g
            b = 0
            [x,y,_] = center_position[i]
            plt.annotate(f'{i}', (x,y), ha='center', color=(r,g,b), weight='bold')
    
    # show plot
    plt.axis('equal')
    plt.show()


## ===============
## Constraints
## ===============
import numpy as np
def set_actives_dep_weights(numUmbrellas, *active_cells, dep_factors = np.logspace(-4, 0, 5)):
    weights_per_cell = np.zeros(numUmbrellas)
    weights_per_cell[active_cells] = 1
    return np.einsum('i, j -> ij', dep_factors, weights_per_cell)

def set_target_height(numUmbrellas, active_cell, target_heights):
    target_height_multiplier = np.ones(numUmbrellas)*-1
    target_height_multiplier[active_cell] = target_heights
    return target_height_multiplier

def percent_to_height(curr_um, thickness, indexes, percents):
    return [((percent/100)*(curr_um.umbrellaHeights[idx]-thickness)+thickness)/thickness
            for percent,idx in zip(percents,indexes)]


import sys
sys.path.append('..')
import umbrella_mesh
def get_arm_stresses(curr_um, stresses):
    stresses = np.array(stresses)
    for sid in range(curr_um.numSegments()):
        seg = curr_um.segment(sid)
        if seg.segmentType() == umbrella_mesh.SegmentType.Plate:
            stresses[sid] = 0
    return stresses.tolist()