import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter # for unique element in surroundings of bumps

import sys
sys.path.append('..')
import configuration

from deployment import deploy
from tools import create_dir_hierarchy

TARGET_PERCENT = 100 # bad habit of global variable, but messy to get it down to `_get_percents`

def deploy_path(input_path, nb_steps, category, surroundings=True, level=1, drop_extrems_at_boundary=False, drop_boudary=False, force_smalest_dep=False, target_percent=100, gif_duration=4, strategies=None, baseline=True, verbose=True):
    
    global TARGET_PERCENT
    TARGET_PERCENT = target_percent
    
    if strategies==None: strategies=get_path_strategies()
    nb_strategies = len(strategies)
        
    io, input_data, target_mesh, curr_um, plate_thickness_scaled, target_height_multiplier = \
        configuration.parse_input(input_path, handleBoundary = False, isHex = False, use_target_surface = False)

    graph = create_graph(input_data['umbrella_connectivity'], curr_um)

    bumps, depressions = find_extrems(graph, surroundings=surroundings, level=level, drop_extrems_at_boundary=drop_extrems_at_boundary, drop_boudary=drop_boudary)
    paths = shortes_paths(graph, bumps, depressions, force_smalest_dep=force_smalest_dep)
    
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
    
def find_extrems(graph, surroundings=True, level=1, drop_extrems_at_boundary=False, drop_boudary=False, degree=3):
    if (surroundings and level<1): raise ValueError("surrounding level should be positive int")
    
    bumps = []
    depression = []
    
    boundary = get_boundary(graph)
    
    graph_copy = graph.copy()
    if drop_boudary: graph_copy.remove_nodes_from(boundary)
    
    for node in graph_copy.nodes:
        is_bump = True
        is_depression = True
        if drop_extrems_at_boundary and graph_copy.degree(node)<degree: continue
        neighbors = [s0 for s in surround_bumps(graph, [node], level=level)[:-1] for s0 in s] if surroundings \
                    else graph_copy.neighbors(node)
        # surround_bumps: -1 is bumps, -2 is 1st level surroundings, -3 is 2nd level... need to be flaten
        
        # if compute 'surround_bumps' with graph_copy, we can have some node with only one link as `bump`
        # and then, the surroundings are not complete... | As it do not alter the graph, it's totally fine
        
        for neighbor in neighbors:
            if drop_boudary and neighbor in boundary : continue
            if graph.nodes[neighbor]['height'] > graph.nodes[node]['height']: is_bump = False
            else : is_depression = False
        if (is_bump): bumps.append(node)
        elif (is_depression): depression.append(node)
    return bumps, depression

def surround_bumps(graph, bumps, level=1, verbose=False, pos=None,
                   min_size=100, max_size=600, colors_default=None, with_labels=False):
    def surround_neigh(graph, curr_surr, prev_surr):
        # first level is quite different
        if len(curr_surr) == 1:
            center = curr_surr[0]
            # 1st degre neighbors (direct link to center)
            neighbors_1 = [n1 for n1 in graph.neighbors(center)]
            # 2nd degre neighbors
            neighbors_2 = []
            for n1 in neighbors_1:
                neighbors_2.extend([n2 for n2 in graph.neighbors(n1) if n2 != center])
            # 3rd degre neighbors
            neighbors_3 = []
            for n2 in neighbors_2:
                for nn2 in graph.neighbors(n2):
                    if nn2 in neighbors_1: continue
                    if nn2 in neighbors_3: continue
                    for nnn2 in graph.neighbors(nn2):
                        if nnn2 == n2: continue
                        if nnn2 in neighbors_2: neighbors_3.append(nn2)
            next_surr = neighbors_1 + neighbors_2 + neighbors_3
            return next_surr, curr_surr
        
        # Level greater than 1
        neighbors = []
        for cs in curr_surr:
            neighbors.extend([n for n in graph.neighbors(cs) if n not in curr_surr and n not in prev_surr])
        next_surr = []
        for n in neighbors:
            next_surr.append(n)
            next_surr.extend([nn for nn in graph.neighbors(n) if nn not in curr_surr
                                                              and nn not in next_surr
                                                              and not graph.nodes.data()[nn]['visited']])
        return next_surr, curr_surr
    
    dic = {}
    for b in bumps:
        dic[b] = [b], []
    shared_graph = graph.copy()
    nx.set_node_attributes(shared_graph, False, 'visited')
    nx.set_node_attributes(shared_graph, {key:True for key in bumps}, 'visited') # set bumps as visited
    
    surr_levels = [bumps]
    for l in range(level):
        surr_level = []
        for key,(curr, prev) in dic.items():
            c, p = surround_neigh(shared_graph, curr, prev)
            surr = [c0 for c0 in c if not shared_graph.nodes.data()[c0]['visited']]
            nx.set_node_attributes(shared_graph, {key:True for key in surr}, 'visited')
            surr_level.extend(surr)
            dic[key] = surr, p
        if surr_level==[] : break # no more medium to propagate the wave
        if (verbose): draw_height_extrems(graph, pos, bumps, surr_level, min_size=min_size, max_size=max_size, colors_default=colors_default, with_labels=with_labels)
        surr_levels.append(surr_level)
    
    return list(reversed(surr_levels)) # bumps are last steps


def shortes_paths(graph, bumps, depressions, force_smalest_dep=False):
    paths = []
    if force_smalest_dep:
        dep_heights = [(idx, data['height']) for (idx, data) in dict(graph.nodes.data()).items() if idx in depressions]
        dep_min = min(dep_heights, key=lambda x: x[1])[0]
        
        shortest = []
        for b in bumps:
            shortest.append(nx.shortest_path(graph,source=b,target=dep_min, weight='weight'))
        lst = min(shortest, key=len)
        paths.append(lst)
        bumps = [b for b in bumps if b!=lst[0]] # from bump to depression -> bump it 1st item
        
    for b in bumps:
        shortest = []
        for d in depressions:
            shortest.append(nx.shortest_path(graph,source=b,target=d, weight='weight'))
        lst = min(shortest, key=len)
        paths.append(lst)
    return paths

def get_boundary(graph, degree=3):
    # return [n for n in graph.nodes if graph.degree(n)<degree]
    return [n[0] for n in graph.degree() if n[1]<degree]

# ======================================================================
# === DRAW ===
# ======================================================================
'''
# light_grey = 211/255
# colors_default = [(light_grey,light_grey,light_grey)] * len(graph) # light_grey
'''
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
    

# ======================================================================
# === Strategies ===
# ======================================================================
def get_path_strategies():
    return *_get_paths_strategies(),    \
           *_get_bumps_strategies(),    \
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


# ----------------
# --- BOUNDARY ---
# ----------------
def _boundary_units(graph, nb_steps, degree=3):
    active_units = get_boundary(graph, degree)
    return active_units,\
           *_get_per_stp(active_units, nb_steps)

def increasing_height(graph):
    increasing_heights = [x[0] for x in sorted(graph.nodes.data(), key=lambda x: x[1]['height'])]
    active_units = [increasing_heights[:i+1] for i in range(len(increasing_heights))] # cumulative
    nb_steps = len(active_units)
    return active_units,\
           *_get_per_stp(active_units, nb_steps)
    
# ======================================================================
# === HELPERS ===
#======================================================================
def _get_percents(pahts):
    return [[TARGET_PERCENT]*len(p) for p in pahts]

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