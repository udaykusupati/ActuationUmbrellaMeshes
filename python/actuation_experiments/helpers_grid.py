import numpy as np

# ======================================================================
# =============================================================== GRID =
# ======================================================================
def set_actives_dep_weights(numUmbrellas, *active_cells,
                            dep_factors = np.logspace(-4, 0, 5)):
    weights_per_cell = np.zeros(numUmbrellas)
    weights_per_cell[active_cells] = 1
    return np.einsum('i, j -> ij', dep_factors, weights_per_cell)

def set_target_height(numUmbrellas, active_cell, target_heights):
    target_height_multiplier = np.ones(numUmbrellas)*-1
    target_height_multiplier[active_cell] = target_heights
    return target_height_multiplier

def percent_to_height(init_height, thickness, indexes, percents):
    return [((1-percent/100)*(init_height[idx]-thickness)+thickness)/thickness
            for percent,idx in zip(percents,indexes)]