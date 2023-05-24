import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter # for unique element in surroundings of bumps

import sys
sys.path.append('..')
import configuration

from deployment import deploy
from tools import create_dir_hierarchy, gif_to_img_duration


def deploy_path(input_path, nb_steps, category, gif_duration=4, strategies=None, baseline=True, verbose=True):
       
    if strategies==None: strategies=get_strategies()
    nb_strategies = len(strategies)
        
    io, input_data, target_mesh, curr_um, plate_thickness_scaled, target_height_multiplier = \
        configuration.parse_input(input_path, handleBoundary = False, isHex = False, use_target_surface = False)

    graph = create_graph(input_data['umbrella_connectivity'], curr_um)

    bumps, depressions = find_extrems(graph)
    surrounds_bumps = surround_bumps(graph, bumps)
    paths = shortes_paths(graph, bumps, depressions)
    
    # check that all strategies are compatible with nb_steps
    for strat  in strategies:
        strat(paths, nb_steps)
    
    
    degree = 3
    rows = 0
    cols = 0
    deployment = 'linear'
    init_heights = curr_um.umbrellaHeights
    
    saved_paths = []
    
    # BASELINE - active all cells:
    if baseline:
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
               active_units, target_percents, init_heights, plate_thickness_scaled,
               gif_duration, deployment, verbose)
        
    # Other Stragegies
    for i, strat in enumerate(strategies):
        if verbose: print("\n\n=====")
        if verbose: print(   f"===== {strat.__name__[1:]} =====")
        if verbose: print(   f"===== {i+1:0>2}/{nb_strategies:0>2}")
        
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
               active_units, target_percents, init_heights, plate_thickness_scaled,
               gif_duration, deployment, verbose)
    
    return saved_paths    
    
def create_graph(connectivity, curr_um):
    graph = nx.Graph()
    
    height_as_node = []
    for i,h in enumerate(curr_um.umbrellaHeights):
        height_as_node.append((i,{'height': h}))
    graph.add_nodes_from(height_as_node)
    
    wighted_edge = []
    for i,j in connectivity:
        wighted_edge.append((i,j,abs(curr_um.umbrellaHeights[i] - curr_um.umbrellaHeights[j])))
    graph.add_weighted_edges_from(wighted_edge)
    
    return graph
    
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

def surround_bumps(graph, bumps):
    'broken if drop_boundary=False'
    graph_copy = graph.copy()
    neighbors = []
    for bump in bumps:
        neighbors.append([n for n in graph_copy.neighbors(bump)])
        graph_copy.remove_node(bump)
    
    surround_neighbors = {}
    degree = 3
    neigh_path_len = 5
    for i,(b,n) in enumerate(zip(bumps,neighbors)):
        paths = []
        for d in range(degree):
            source=n[d]
            target=n[(d+1)%degree]
            # are they boundary ? (the shortest path thus goes through the 3rd neighbor)
            if len(list(graph_copy.neighbors(source)))==1  and\
               len(list(graph_copy.neighbors(target)))==1: continue
            paths.append(nx.shortest_path(graph_copy, source=source, target=target, weight=None))
        
        # get rid of duplicate units
        init_path = 0
        if len(paths)==degree : init_path = 1
        surround_neighbors[b] = paths[0][init_path:]
        for path in paths[1:]:
            surround_neighbors[b].extend(path[1:])
    return surround_neighbors

def shortes_paths(graph, bumps, depressions):
    paths = []
    for b in bumps:
        shortest = []
        for d in depressions:
            shortest.append(nx.shortest_path(graph,source=b,target=d, weight='weight'))
        lst = min(shortest, key=len)
        paths.append(lst)
    return paths

# ======================================================================
# === DRAW ===
# ======================================================================
def draw_height(graph, pos, min_size=100, max_size=600, with_labels=False):
    scaled_height = _scaled_height(graph, min_size=min_size, max_size=max_size)
    nx.draw(graph, pos=pos, with_labels=with_labels, node_size=scaled_height)
    plt.show()
    
def draw_height_extrems(graph, pos, bumps, depressions, min_size=100, max_size=600, colors_default=None, with_labels=False):
    colors = colors_default.copy() if not colors_default==None else ['#1f78b4']*len(graph)# '#1f78b4' is nx default value
    for bump in bumps:
        colors[bump] = (1,0,0)
    for dep in depressions:
        colors[dep] = (0,1,0)
    _draw_colored(graph, pos, colors, min_size=min_size, max_size=max_size, with_labels=with_labels)
    
def draw_height_path(graph, pos, paths, min_size=100, max_size=600, colors_default=None, with_labels=False):
    colors = colors_default.copy() if not colors_default==None else ['#1f78b4']*len(graph)# '#1f78b4' is nx default value
    for path in paths:
        path_length = len(path)-1
        for i, node in enumerate(path): # path is from bump to depression
            g = i/path_length
            r = 1-g
            b = 0
            colors[node] = (r,g,b)
    _draw_colored(graph, pos, colors, min_size=min_size, max_size=max_size, with_labels=with_labels)
    
def draw_height_surrounds(graph, pos, surrounds, min_size=100, max_size=600, colors_default=None, with_labels=False):
    colors = colors_default.copy() if not colors_default==None else ['#1f78b4']*len(graph)# '#1f78b4' is nx default value
    
    for s in surrounds.items():
        colors[s[0]] = (1,0,0) # bumps in red
        for neighbor in s[1]:
            colors[neighbor] = (0,1,0) # bumps in red
    _draw_colored(graph, pos, colors, min_size=min_size, max_size=max_size, with_labels=with_labels)
    

# ======================================================================
# === Strategies ===
# ======================================================================
def get_strategies():
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


# ======================================================================
# === HELPERS ===
#======================================================================
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

def _scaled_height(graph, min_size=100, max_size=600):
    heights = [attr[1]['height'] for attr in graph.nodes.data()]
    min_h, max_h = min(heights), max(heights)
    factor = (max_size-min_size)/(max_h-min_h)
    return [(h-min_h)*factor+min_size for h in heights]

def _draw_colored(graph, pos, colors, min_size=100, max_size=600, with_labels=False):
    scaled_height = _scaled_height(graph, min_size=min_size, max_size=max_size)
    nx.draw(graph, pos=pos, with_labels=with_labels, node_size=scaled_height, node_color=colors)
    plt.show()