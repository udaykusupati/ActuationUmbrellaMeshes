import sys
import os
sys.path.append('..')
import umbrella_mesh

import numpy as np

# ======================================================================
# ============================================================== TOOLS =
# ======================================================================
def create_dir(name):
    if not os.path.exists(name): os.makedirs(name)
    else: raise ValueError(f'folder {name} already exists')

def get_max_stress_matrix(curr_um,
                          stress_type='maxBending'):
    '''return the adjacency max stress matrix'''
    matrix = np.zeros((curr_um.numUmbrellas(), curr_um.numUmbrellas()))
    stress_all = get_stresses(curr_um, stress_type)
    for seg_id in range(curr_um.numSegments()):
        seg = curr_um.segment(seg_id)
        if seg.segmentType()==umbrella_mesh.SegmentType.Arm:
            sjoint = seg.startJoint
            ejoint = seg.endJoint
            stresses = stress_all[seg_id]
            if curr_um.joint(sjoint).jointPosType()==umbrella_mesh.JointPosType.Arm:
                i,j = curr_um.joint(ejoint).umbrellaID()
            if curr_um.joint(ejoint).jointPosType()==umbrella_mesh.JointPosType.Arm:
                i,j = curr_um.joint(sjoint).umbrellaID()
            matrix[i,j] = matrix[j,i] = max(matrix[i,j], max(stresses))
    return matrix

def get_stresses(curr_um,
                 stress_type):
    if stress_type=='VonMises':
        stresses = curr_um.maxVonMisesStresses()
    elif stress_type=='maxBending':
        stresses = curr_um.maxBendingStresses()
    elif stress_type=='minBending':
        stresses = curr_um.minBendingStresses()
    elif stress_type=='Twisting':
        stresses = curr_um.twistingStresses()
    elif stress_type=='Stretching':
        stresses = curr_um.stretchingStresses()
    else:
        raise ValueError(f'the required stress type <{stress_type}> do not correspond to any available stress')
    return stresses

def get_stresses_types():
    return ['VonMises',
            'maxBending', 
            'minBending',
            'Twisting',
            'Stretching']

def get_smatrix_min_max(curr_um, stress_type,
                         zero_as_extrem=False):
    matrix = get_max_stress_matrix(curr_um, stress_type)
    if not zero_as_extrem:
        tmp = matrix[matrix != 0]
        return matrix, tmp.min(), tmp.max()
    return matrix, matrix.min(), matrix.max()


def set_plate_stress_null(curr_um,
                           stresses):
    stresses = np.array(stresses)
    for sid in range(curr_um.numSegments()):
        seg = curr_um.segment(sid)
        if seg.segmentType() == umbrella_mesh.SegmentType.Plate:
            stresses[sid] = 0
    return stresses.tolist()