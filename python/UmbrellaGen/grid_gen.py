import sys
sys.path.append('UmbrellaGen/')
import umbrella_gen
import utils
import numpy as np
import math

def genUmbrellaWithHeights(degree = 3, rows = 5, cols = 5, height_scales = None, minHeight = 64, useOverHang = False, armPlateEdgeAxisOffset = 10, armJointAxisOffset = 10, minOverHang = 5.0, plate_thickness = 3.0, edge_length = 30, base_mesh = None, scaling = True):
    edgeLength = edge_length
    t_mesh = None
    uv_mesh = None
    # max_overhang = edgeLength
    if height_scales is not None:
        max_overhang = np.max(np.array(height_scales))*minHeight - minHeight
        fabHeight = minHeight + max_overhang
    jsonPath = '../UmbrellaGen/grid_dump.json.gz'

    if base_mesh is not None:
        # base_mesh has to be equilateral
        v_base, f_base = base_mesh
        v_out, f_out = v_base.copy(), f_base.copy()
        V, F = np.array(v_out), np.array(f_out)
        c_out = V[F].mean(axis=1)

        assert np.unique(np.round(np.linalg.norm(np.expand_dims(c_out, 1) - V[F], axis = -1), decimals = 7)).size == 1, "Base mesh has to be equilateral"
        c_out = c_out.tolist()
        x_out = [[] for _ in F]
        for fid1, f1 in enumerate(F):
            for fid2, f2 in enumerate(F):
                if len(np.intersect1d(f1, f2)) == 2:
                    if fid2 not in x_out[fid1]: x_out[fid1].append(fid2)
                    if fid1 not in x_out[fid2]: x_out[fid2].append(fid1)
        i_out = None
    elif degree == 3:
        v_out, f_out, c_out, i_out, x_out = utils.generate_regular_tri_grid(rows + 1, cols + 1, l = edgeLength)
    elif degree == 4:
        v_out, f_out, c_out, i_out, x_out = utils.generate_regular_sqr_grid(rows + 1, cols + 1, l = edgeLength)
    elif degree == 6:
        v_out, f_out, c_out, i_out, x_out = utils.generate_regular_tri_grid(rows, cols, l = edgeLength * utils.sqrt(3)) # Scale up the dual mesh so that the hex has same side length as other cells
    elif degree == 99: ## SPECIAL HARD CODE FOR ONE TOPOLOGY
        degree = 3
        v_out, f_out, c_out, i_out, x_out = utils.generate_regular_tri_grid(rows + 1, cols + 1, l = edgeLength)
        invalid_tris = [0,1,2, 8,9,10,11,19,20,30,40,41,49]
        utils.noneList(f_out, invalid_tris)
        utils.noneList(c_out, invalid_tris)
        utils.noneList(i_out, invalid_tris)
        utils.noneList(x_out, invalid_tris)
        c_out, _ = utils.noneListR(c_out)
        i_out, _ = utils.noneListR(i_out)
        x_out, _ = utils.noneListR(x_out)
        f_out, idxD = utils.noneListR(f_out)

        temp = []
        for xi in x_out:
            tempL = [idxD[xii] for xii in xi if xii in idxD]
            temp.append(tempL)
        x_out = temp

        
    else:
        assert 0 
    v_ret, f_ret, x_ret, c_ret = v_out.copy(), f_out.copy(), x_out.copy(), c_out.copy()
    
    scaleLength = height_scales
    overhangs = None
    if height_scales is not None:
        overhangs = [0.0] * len(f_out)
        assert len(height_scales) == len(f_out) or len(height_scales) == len(v_out)
        height_scales = np.array(height_scales)
        if len(height_scales[height_scales < 1]) > 0: assert 0, "Height must be atleast min Height"
        # if len(height_scales[height_scales > 1 + max_overhang/minHeight]) > 0: assert 0, ("Max Height achievable is min Height + max_overhang, so max height scale is: ", 1 + max_overhang/minHeight)
        overhangs = (fabHeight - height_scales * minHeight + minOverHang).tolist() # Fix min overhang. Currently set to 5.0
        # overhangs = [10.0] * len(f_out)
        # print(overhangs)
        heights = height_scales * minHeight
    else:
        scaleLength = [1]*len(f_out)
    if not useOverHang:
        overhangs = None
        heights = None
    umbrella_gen.genPattern(edgeLength, t_mesh, uv_mesh, i_out, v_out, f_out, c_out, x_out, scaleLength, jsonPath, marginLength = 0.0, armPlateEdgeAxisOffset = armPlateEdgeAxisOffset, armJointAxisOffset = armJointAxisOffset, asymmetryOffset = 0, width = 5, thickness = plate_thickness, targetSpacingFactor = 5, minHeight = minHeight, select_umbrella = False, handlePivots = False, min_coerce_dist = 1e-8, degree = degree, overhang=overhangs, heights = heights, scaling = scaling)
    # return v_ret, f_ret, x_ret, c_ret
    


def scaleFactor(m, h, s):
    return 1 + math.sqrt(3*h**2 - 3*s**2)/m
def spacingFromCF(m, h, c):
    return math.sqrt(h**2 - ((m*(c - 1))**2)/3)

def maximize_deployment(old_heights, quantizedHeights, deployment_spacings, plateEdgeLength):
    new_deployment_spacings = np.zeros_like(deployment_spacings)
    new_heights = np.zeros_like(old_heights)
    print(old_heights)
    print(deployment_spacings)
    for id, old_height in enumerate(old_heights):
        new_height_id = ((quantizedHeights - deployment_spacings[id])/quantizedHeights).argmax()
        new_height = quantizedHeights[new_height_id]
        conformal_factor = scaleFactor(plateEdgeLength, old_height, deployment_spacings[id])
        new_deployment_spacings[id] = spacingFromCF(plateEdgeLength, new_height, conformal_factor)
        new_heights[id] = new_height
    return new_heights, new_deployment_spacings
        


