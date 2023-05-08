import networkx as nx

import sys
sys.path.append('..')
import configuration

from deployment import deploy


# deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
#             active_cells, target_percents, init_heights, plate_thickness, deployment,
#             verbose=False)

'''
- generate mesh
- generate graph
- get shortest path
- deploy
'''

def deploy_path(input_path):
    
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
    
    # deploy(path, folder_name, input_data, curr_um, degree, rows, cols, steps,
    #         active_cells, target_percents, init_heights, plate_thickness_scaled, deployment,
    #         verbose)
    
    return paths
    
    
    
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

    