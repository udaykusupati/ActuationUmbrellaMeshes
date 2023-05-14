import networkx as nx

import sys
sys.path.append('..')
import configuration

from deployment import deploy
from tools import create_dir_hierarchy

# deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
#             active_cells, target_percents, init_heights, plate_thickness, deployment,
#             verbose=False)

'''
- generate mesh
- generate graph
- get shortest path
- deploy
'''

def deploy_path(input_path, nb_steps, category, verbose=True):
    
    io, input_data, target_mesh, curr_um, plate_thickness_scaled, target_height_multiplier = \
        configuration.parse_input(input_path, handleBoundary = False, isHex = False, use_target_surface = False)

    graph = nx.Graph()
    
    height_as_node = []
    for i,h in enumerate(curr_um.umbrellaHeights):
        height_as_node.append((i,{'height': h}))
    graph.add_nodes_from(height_as_node)
    
    wighted_edge = []
    for i,j in input_data['umbrella_connectivity']:
        wighted_edge.append((i,j,abs(curr_um.umbrellaHeights[i] - curr_um.umbrellaHeights[j])))
    graph.add_weighted_edges_from(wighted_edge)

    paths = shortes_paths(graph, *find_extrems(graph))
    
    # aps : active_units, target_percents, steps
    aps_strategies = [s(paths, nb_steps) for s in _get_strategies()]
        
    degree = 3
    rows = 0
    cols = 0
    deployment = 'linear'
    init_heights = curr_um.umbrellaHeights
    
    saved_paths = []
    
    # BASELINE - active all cells:
    strat = _units_all_at_once
    if verbose: print(f"\n\n===== {strat.__name__[1:]} =====\n")
    active_units, target_percents, steps = strat(paths, curr_um.numUmbrellas(), nb_steps)
    folder_name , path = create_dir_hierarchy(category,
                                              degree,
                                              rows,
                                              cols,
                                              deployment,
                                              strat.__name__[1:])
    saved_paths.append(path)

    deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
            active_units, target_percents, init_heights, plate_thickness_scaled, deployment,
            verbose)
        
    strats = _get_strategies()
    for i, strat in enumerate(strats):
        if verbose: print("\n\n=====")
        if verbose: print(   f"===== {strat.__name__[1:]} =====")
        if verbose: print(   f"===== {i+1:0>2}/{len(strats):0>2}")
        
        # re-generate curr_um for undeployed plot...
        io, input_data, target_mesh, curr_um, plate_thickness_scaled, target_height_multiplier = \
        configuration.parse_input(input_path, handleBoundary = False, isHex = False, use_target_surface = False)
        
        active_units, target_percents, steps = strat(paths, nb_steps)
        folder_name , path = create_dir_hierarchy(category,
                                                  degree,
                                                  rows,
                                                  cols,
                                                  deployment,
                                                  strat.__name__[1:])
        saved_paths.append(path)
        
        deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
                active_units, target_percents, init_heights, plate_thickness_scaled, deployment,
                verbose)
    
    return saved_paths
    
    
    
def find_extrems(graph, drop_boudary=True, degree=3):
    bumps = []
    depression = []
    for node in graph.nodes:
        is_bump = True
        is_depression = True
        if drop_boudary and len(list(graph.neighbors(node)))<degree: continue
        for neighbor in graph.neighbors(node):
            if graph.nodes[neighbor]['height'] > graph.nodes[node]['height']: is_bump = False
            else : is_depression = False
        if (is_bump): bumps.append(node)
        elif (is_depression): depression.append(node)
    return bumps, depression

def shortes_paths(graph, bumps, depressions):
    paths = []
    for b in bumps:
        shortest = []
        for d in depressions:
            shortest.append(nx.shortest_path(graph,source=b,target=d))
        lst = min(shortest, key=len)
        paths.append(lst)
    return paths


# === Strategies ===

def _get_strategies():
    return *_get_paths_strategies(),\
           *_get_bumps_strategies(),\
           *_get_actuators_strategies()


# -----------
# --- ALL ---
# -----------
def _units_all_at_once(paths, nb_units, nb_steps):
    active_units = [list(range(nb_units))]
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

# -------------
# --- PATHS ---
# -------------
def _paths_all_at_once(paths, nb_steps):
    active_units = [[unit for units in paths for unit in units]]
    return active_units,\
           *_get_per_stp(active_units, nb_steps)
    
def _paths_one_after_one(paths, nb_steps):
    active_units = paths.copy()
    return active_units,\
           *_get_per_stp(active_units, nb_steps)
    
def _paths_cum(paths, nb_steps):
    active_units = []
    temp = []
    for p in paths:
        temp.extend(p)
        active_units.append(temp.copy())
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

# -------------------------
def _get_paths_strategies():
    return _paths_all_at_once,\
           _paths_one_after_one,\
           _paths_cum


# -------------
# --- BUMPS ---
# -------------
def _bumps_all_at_once(paths, nb_steps):
    active_units = [[path[0] for path in paths]]
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

def _bumps_one_after_one(paths, nb_steps):
    active_units = [[path[0]] for path in paths]
    return active_units,\
           *_get_per_stp(active_units, nb_steps)
    
def _bumps_cum(paths, nb_steps):
    # cummulative bumps
    active_units = []
    temp = []
    for p in paths:
        temp.extend([p[0]])
        active_units.append(temp.copy())
        
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

# -------------------------
def _get_bumps_strategies():
    return _bumps_all_at_once,\
           _bumps_one_after_one,\
           _bumps_cum


# -----------------
# --- ACTUATORS ---
# -----------------
def _actuators_path_one_after_one(paths, nb_steps):
    rpahts = _reverse_path(paths)
    
    active_units = []
    for path in rpahts:
        for a,b in zip(path[:-1], path[1:]):
            active_units.append([a,b])
            
    return active_units,\
           *_get_per_stp(active_units, nb_steps)
            
def _actuators_two_per_path(paths, nb_steps):
    rpahts = _reverse_path(paths)
    
    active_units = []
    max_path_length = len(max(paths, key=len))-1
    for i in range(max_path_length):
        temp = []
        for path in rpahts:
            idx = min(len(path)-2, i)
            temp.extend(path[idx:idx+2])
        active_units.append(temp)
        
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

def _actuators_two_per_path_end(paths, nb_steps):
    # all paths 'ends' together
    rpahts = _reverse_path(paths)

    active_units = []
    max_path_length = len(max(paths, key=len))-1
    for i in range(max_path_length):
        temp = []
        for path in rpahts:
            idx = max(i-(max_path_length-len(path))-1, 0)
            temp.extend(path[idx:idx+2])
        active_units.append(temp)
        
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

# -----------------------------
def _get_actuators_strategies():
    return _actuators_path_one_after_one,\
           _actuators_two_per_path,\
           _actuators_two_per_path_end


# ---------------
# --- HELPERS ---
# ---------------
def _get_percents(pahts):
    return [[100]*len(p) for p in pahts]

def _get_steps(pahts, nb_steps):
    nb_phase = len(pahts)
    if nb_phase > nb_steps: raise ValueError(f"required number of phase [{nb_phase}] is greater than the requested number of steps [{nb_steps}]")
    steps = []
    for i in range(nb_phase):
        steps.append(nb_steps//nb_phase)
        
    for i in range(nb_steps%nb_phase):
        steps[-1-i] += 1 # adding more steps to the lasts phases
    return steps

def _get_per_stp(pahts, nb_steps):
    return _get_percents(pahts),\
           _get_steps(pahts, nb_steps)

def _reverse_path(paths):
    reversed_paths = []
    for path in paths:
        reversed_paths.append([x for x in reversed(path)])
    return reversed_paths