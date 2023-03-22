from dis import dis
from hashlib import new
from ntpath import join
import json, copy, gzip
from utils import *

def genLinkageVerts(zHeight, zHeight_nb, uid, nbr_map, c, c_nb, edgeLength, marginLength, armPlateEdgeAxisOffset, armJointAxisOffset, asymmetryOffset, degree = 3, overhang = None):

    # Remember: pivot offset handling would move the TA joints of neighboring plates inward so that a spacing of "thickness" is made between them. Similarly for the z direction. Here the temp plate is constructed which will be made smaller later to truePlateLength = tempPlateLength - 2*thickness*sqrt(3)
    if degree == 3:
        tempPlateLength = edgeLength - 2*sqrt(3.)*marginLength
    elif degree == 4:
        tempPlateLength = edgeLength - 2*marginLength
    elif degree == 6:
        tempPlateLength = edgeLength - 2*marginLength/sqrt(3.)

    halfZheight = zHeight/2.
    if degree == 3:
        dist  = sqrt(3.)/6.*tempPlateLength
    elif degree == 4:
        dist  = tempPlateLength/2.
    elif degree == 6:
        dist  = sqrt(3.)/2.*tempPlateLength

    if len(c_nb) > 0:
        x_axis = normalizeVec(vecSub(c_nb[0], c))
    else:
        x_axis = [1.0,0,0]
    assert x_axis[2] == 0
    assert c[2] == 0
    y_axis = [-x_axis[1], x_axis[0], 0]
    

    v = []

    # Origin Top
    v.append({'type':'OT', 'uid' : uid, 'edges' : [], 'coord': [c[0],c[1],halfZheight]})
    # Origin Bot
    v.append({'type':'OB', 'uid' : uid, 'edges' : [], 'coord': [c[0],c[1],-halfZheight + asymmetryOffset]})
    
    # Top Plate
    v_top, v_top_normals = [], []
    coord = vecAdd(vecMult(x_axis, dist), vecMult(y_axis, armPlateEdgeAxisOffset))
    v_top.append({'type': 'PT', 'edges': [], 'coord': [coord[0], coord[1], halfZheight]})
    v_top_normals.append({'type': 'PTN', 'edges': [], 'coord': [-y_axis[0], -y_axis[1], 0]})
    v_top = [ptRotPoly(vi, degree = degree) for vi in v_top]
    v_top_normals = [ptRotPoly(vi, degree = degree) for vi in v_top_normals]
    v_top, v_top_normals = flatten(v_top), flatten(v_top_normals)
    v_top_check = [{'type': vi['type'], 'edges': vi['edges'], 'coord': rnd(vecAdd(vi['coord'], c))} for vi in v_top]
    v_top = [{'type': vi['type'], 'edges': vi['edges'], 'coord': vecAdd(vi['coord'], c), 'ghost_normal': vin['coord']} for vi, vin in zip(v_top, v_top_normals)]
    for top_id in range(degree):
        v_top[top_id]['uid'] = [uid, -1 - top_id]
    v += v_top

    top_link_ids = []
    for nb_id in range(len(zHeight_nb)):
        nb_x_axis = normalizeVec(vecSub(c_nb[nb_id], c))
        nb_y_axis = [-nb_x_axis[1], nb_x_axis[0], 0]
        coord = vecAdd(vecMult(nb_x_axis, dist), vecMult(nb_y_axis, armPlateEdgeAxisOffset)) 
        v_temp = {'type': 'PT', 'edges': [], 'coord': rnd([c[0] + coord[0], c[1] + coord[1], halfZheight])}
        assert v_temp in v_top_check
        if v_temp in v_top_check:
            v_top_id = v_top_check.index(v_temp)
            top_link_ids.append(v_top_id + 2)
    for nb_id, top_id in enumerate(top_link_ids):
        v[top_id]['uid'] = [uid, nbr_map[uid][nb_id]]

    # Bot Plate


    v_bot, v_bot_normals = [], []
    coord = vecAdd(vecMult(x_axis, dist), vecMult(y_axis, -armPlateEdgeAxisOffset))
    v_bot.append({'type': 'PB', 'edges': [], 'coord': [coord[0], coord[1], -halfZheight + asymmetryOffset]})
    v_bot_normals.append({'type': 'PBN', 'edges': [], 'coord': [y_axis[0], y_axis[1], 0]})
    v_bot = [ptRotPoly(vi, degree = degree) for vi in v_bot]
    v_bot_normals = [ptRotPoly(vi, degree = degree) for vi in v_bot_normals]
    v_bot, v_bot_normals = flatten(v_bot), flatten(v_bot_normals)
    v_bot_check = [{'type': vi['type'], 'edges': vi['edges'], 'coord': rnd(vecAdd(vi['coord'], c))} for vi in v_bot]
    v_bot = [{'type': vi['type'], 'edges': vi['edges'], 'coord': vecAdd(vi['coord'], c), 'ghost_normal': vin['coord']} for vi, vin in zip(v_bot, v_bot_normals)]
    for bot_id in range(degree):
        v_bot[bot_id]['uid'] = [uid, -1 - bot_id]
    v += v_bot
    for nb_id, top_id in enumerate(top_link_ids):
        v[top_id + degree]['uid'] = [uid, nbr_map[uid][nb_id]]

    
    # Arm-Arm Joint
    
    for nb_id in range(len(zHeight_nb)):
        zHeight_nbr = zHeight_nb[nb_id]
        nb_x_axis = normalizeVec(vecSub(c_nb[nb_id], c))
        nb_y_axis = [-nb_x_axis[1], nb_x_axis[0], 0]
        coord = vecAdd(vecMult(nb_x_axis, dist + 2*marginLength*zHeight/(zHeight + zHeight_nbr)), vecMult(nb_y_axis, armJointAxisOffset)) 
        v.append({'type': 'AA', 'edges': [], 'coord': [c[0] + coord[0], c[1] + coord[1], 0]})
        if  True: #uid < nbr_map[uid][nb_id]:
            v[-1]['uid'] = [uid, nbr_map[uid][nb_id]]
        else:
            v[-1]['uid'] = [nbr_map[uid][nb_id], uid]
    # Neighbor-Plate Joint Bot
    for nb_id in range(len(zHeight_nb)):
        zHeight_nbr = zHeight_nb[nb_id]

        nb_x_axis = normalizeVec(vecSub(c_nb[nb_id], c))
        nb_y_axis = [-nb_x_axis[1], nb_x_axis[0], 0]
        coord = vecAdd(vecMult(nb_x_axis, dist + 2*marginLength), vecMult(nb_y_axis, armPlateEdgeAxisOffset)) 
        v.append({'type': 'PB', 'uid': [nbr_map[uid][nb_id], uid], 'edges': [], 'coord': [c[0] + coord[0], c[1] + coord[1], -zHeight_nbr/2 + asymmetryOffset], 'ghost_normal': [-nb_y_axis[0], -nb_y_axis[1], 0]})

     # Top and Bot Overhang
    if overhang is not None:
        ## Straight arms 
        assert armJointAxisOffset == armPlateEdgeAxisOffset
        for nb_id in range(len(zHeight_nb)):
            zHeight_nbr = zHeight_nb[nb_id]
            nb_x_axis = normalizeVec(vecSub(c_nb[nb_id], c))
            nb_y_axis = [-nb_x_axis[1], nb_x_axis[0], 0]
            overhang_height = overhang[uid]/2 * (zHeight + zHeight_nbr)/sqrt(((zHeight + zHeight_nbr))**2 + (2 * marginLength)**2)
            # print(uid, overhang[uid]/2, overhang_height)
            coord = vecAdd(vecMult(nb_x_axis, dist - overhang[uid] * (2 * marginLength) / sqrt(((zHeight + zHeight_nbr))**2 + (2 * marginLength)**2)), vecMult(nb_y_axis, armPlateEdgeAxisOffset)) 
            # Top Overhang (TH)
            v.append({'type': 'TH', 'uid': [uid, nbr_map[uid][nb_id] ], 'edges': [], 'coord': [c[0] + coord[0], c[1] + coord[1], halfZheight + overhang_height], 'ghost_normal': [-nb_y_axis[0], -nb_y_axis[1], 0]})
        
        for nb_id in range(len(zHeight_nb)):
            zHeight_nbr = zHeight_nb[nb_id]
            nb_x_axis = normalizeVec(vecSub(c_nb[nb_id], c))
            nb_y_axis = [-nb_x_axis[1], nb_x_axis[0], 0]
            overhang_height = overhang[nbr_map[uid][nb_id]]/2 * (zHeight + zHeight_nbr)/sqrt(((zHeight + zHeight_nbr))**2 + (2 * marginLength)**2)
            # Bot Overhang (BH)
            coord = vecAdd(vecMult(nb_x_axis, dist +  (2 * marginLength) *(1 + overhang[nbr_map[uid][nb_id]] / sqrt(((zHeight + zHeight_nbr))**2 + (2 * marginLength)**2))), vecMult(nb_y_axis, armPlateEdgeAxisOffset)) 
            v.append({'type': 'BH', 'uid': [nbr_map[uid][nb_id], uid], 'edges': [], 'coord': [c[0] + coord[0], c[1] + coord[1], -zHeight_nbr/2 + asymmetryOffset -overhang_height], 'ghost_normal': [-nb_y_axis[0], -nb_y_axis[1], 0]})
    
    return v, top_link_ids
    
def rnd(vec):
    return [round(v, 6) for v in vec]
def r(d):
    return d#round(d,8)
def genLinkage(vL, topL, degree = 3, overhang = None):
    rsEdge = []
    edgeList = []
    uid = 0
    for v, top_link_ids in zip(vL, topL):
        # len(v) = 2 * num_nbr # AA and PB (neighbor)
        #        + 2 * degree # PT and PB
        #        + 2  # OT and OB
        #        + 2 * num_nbr # TH and BH

        num_nbr = (len(v) - (2 + degree * 2))//(2 + 2 * (overhang is not None))
        assert (len(v) - (2 + degree * 2))%(2 + 2 * (overhang is not None)) == 0
        conList = []
        for eid in range(degree):
            conList += [{'type': 'PT', 'uid': uid, 'con':[ 0, 2 + eid]}] # top poly
            conList += [{'type': 'PB', 'uid': uid, 'con':[ 1, 2 + degree + eid]}] # bot poly
        
        for nb_id in range(num_nbr):
            conList += [{'type': 'TA', 'uid': uid, 'con': [top_link_ids[nb_id], 2 + degree * 2 + nb_id]}]
            conList += [{'type': 'ANB', 'uid': uid, 'con': [2 + degree * 2 + nb_id, 2 + degree * 2 + num_nbr + nb_id]}]
            if overhang is not None:
                conList += [{'type': 'HT', 'uid': uid, 'con':[top_link_ids[nb_id], 2 + degree * 2 + num_nbr*2 + nb_id]}]
                conList += [{'type': 'NBH', 'uid': uid, 'con':[2 + degree * 2 + num_nbr + nb_id, 2 + degree * 2 + num_nbr*3 + nb_id]}]
        edgeList.append(conList)
        uid += 1
    return [rsEdge,edgeList]
    
def newPt(vGlb, vGlb_topo, vi):
    vi_topo = {'type': vi['type'], 'uid': vi['uid']}
    if vi_topo not in vGlb_topo:
        vGlb.append(copy.deepcopy(vi))
        vGlb_topo.append(vi_topo)
        idx = len(vGlb) - 1
    else:
        idx = vGlb_topo.index(vi_topo)
    return [vGlb, vGlb_topo, idx]


def collect_all(vList, edgeList):
    vGlb = []
    vGlb_aux = []
    edgeGlb = []
    uid = 0
    for ei, vi in zip(edgeList, vList):
        for eii in ei:
            idx = [-1,-1]
            [vGlb, vGlb_aux, idx[0]] = newPt(vGlb, vGlb_aux, vi[eii['con'][0]]) # index 1
            [vGlb, vGlb_aux, idx[1]] = newPt(vGlb, vGlb_aux, vi[eii['con'][1]]) # index 2
            edgeGlb.append({'type' : eii['type'], 'con': idx, 'uid': eii['uid']})
        uid+=1
    nbr_map = [[] for _ in range(uid)]
    for vid in range(len(vGlb)):
        if vGlb[vid]['type'] == 'PB' or vGlb[vid]['type'] == 'PT' or vGlb[vid]['type'] == 'AA':
            vGlb[vid]['uid'] = [-1, -1] # First one would be the umbrella having this vertex as plate. Second one connecting this to its arm.
    for eid, edge in enumerate(edgeGlb):
        for con_id in [0,1]:
            vGlb[edge['con'][con_id]]['edges'].append(eid)
            if vGlb[edge['con'][con_id]]['type'] == "PB":
                if edge['type'] == 'PB':
                    vGlb[edge['con'][con_id]]['uid'][0] = edgeGlb[eid]['uid']
                elif edge['type'] == 'ANB':
                    vGlb[edge['con'][con_id]]['uid'][1] = edgeGlb[eid]['uid']
                else: assert edge['type'] == 'NBH'
            elif vGlb[edge['con'][con_id]]['type'] == "AA":
                if edge['type'] == 'TA':
                    vGlb[edge['con'][con_id]]['uid'][0] = edgeGlb[eid]['uid']
                # Other end uid will be filled in later once all PBs have both uids
            elif vGlb[edge['con'][con_id]]['type'] == "PT":
                if edge['type'] == 'PT':
                    vGlb[edge['con'][con_id]]['uid'][0] = edgeGlb[eid]['uid']
                # Other end uid will be filled in later once all AAs have both uids
            elif vGlb[edge['con'][con_id]]['type'] == "TH":
                assert edge['type'] == 'HT'
                vGlb[edge['con'][con_id]]['uid'] = edgeGlb[eid]['uid']
            elif vGlb[edge['con'][con_id]]['type'] == "BH":
                assert edge['type'] == 'NBH', edge['type']
                # Will be filled once all PBs have both uids
            else: # OT OB 
                vGlb[edge['con'][con_id]]['uid'] = edgeGlb[eid]['uid']
    
    for vid, v in enumerate(vGlb):
        if v['type'] == 'PB':
            for eid in v['edges']:
                if edgeGlb[eid]['type'] == 'NBH':
                    assert vGlb[edgeGlb[eid]['con'][1]]['type'] == 'BH'
                    vGlb[edgeGlb[eid]['con'][1]]['uid'] = v['uid'][0]
                if edgeGlb[eid]['type'] == 'ANB':
                    assert vGlb[edgeGlb[eid]['con'][0]]['type'] == 'AA'
                    vGlb[edgeGlb[eid]['con'][0]]['uid'][1] = v['uid'][0]

    for vid, v in enumerate(vGlb):
        if v['type'] == 'AA':
            assert -1 not in v['uid']
            for eid in v['edges']:
                if edgeGlb[eid]['type'] == 'TA':
                    assert vGlb[edgeGlb[eid]['con'][0]]['type'] == 'PT'
                    vGlb[edgeGlb[eid]['con'][0]]['uid'][1] = v['uid'][1]

    for vid in range(len(vGlb)):
        if vGlb[vid]['type'] == 'PB':
            if vGlb[vid]['uid'][0] != -1 and vGlb[vid]['uid'][1] != -1:
                if vGlb[vid]['uid'][1] not in nbr_map[vGlb[vid]['uid'][0]]:
                    nbr_map[vGlb[vid]['uid'][0]].append(vGlb[vid]['uid'][1])
                if vGlb[vid]['uid'][0] not in nbr_map[vGlb[vid]['uid'][1]]:
                    nbr_map[vGlb[vid]['uid'][1]].append(vGlb[vid]['uid'][0])
    return vGlb, edgeGlb, nbr_map

def get_edge_vec(vGlb, edgeGlb, eid, oid):
    # oid : id of origin point of returned edge_vec
    edge = edgeGlb[eid]
    if edge['con'][0] == oid:
        p1_id, p2_id = oid, edge['con'][1]
    elif edge['con'][1] == oid:
        p1_id, p2_id = oid, edge['con'][0]
    else:
        assert 0
    p1 = vGlb[p1_id]['coord']
    p2 = vGlb[p2_id]['coord']
    return vecSub(p2, p1)

def collect_ghost_data(vGlb, edgeGlb, degree = 3):
    up = [0.0, 0.0, 1.0]
    down = [0.0, 0.0, -1.0]
    x_axis = [1.0, 0.0, 0.0]
    alpha_tol = 1e-6
    for eid, edge in enumerate(edgeGlb):
        if edge['type'] == 'PT':
            edge['segment_normal'] = up
        if edge['type'] == 'PB':
            edge['segment_normal'] = down
        edgeGlb[eid] = edge
    for vid, v in enumerate(vGlb):
        if v['type'][0] == 'O': # Rigid triangle center joint
            v['alpha'] = 0
            v['is_rigid'] = True
            v['ghost_bisector'] = x_axis
            if v['type'][1] == 'T':
                v['ghost_normal'] = up
            elif v['type'][1] == 'B':
                v['ghost_normal'] = down
            else: assert 0

            v['A_segments'] = v['edges']
            v['B_segments'] = []
            v['midpoint_offsets_A'] = []
            v['midpoint_offsets_B'] = []
            for seg in v['A_segments']:
                v['midpoint_offsets_A'].append([0, 0, 0])
            assert len(v['A_segments']) == degree and len(v['B_segments']) == 0
        
        if v['type'][0] == 'P': # Plate to Arm Joints 

            # Ghost Normal already collected
            # Mind that edge PT isn't normal to the cell edge
            
            v['A_segments'] = []
            v['B_segments'] = []
            for eid in v['edges']:
                if edgeGlb[eid]['type'] == 'TA' or edgeGlb[eid]['type'] == 'ANB' or edgeGlb[eid]['type'] == 'HT' or edgeGlb[eid]['type'] == 'NBH':
                    v['B_segments'].append(eid)
                if edgeGlb[eid]['type'][0] == 'P':
                    v['A_segments'].append(eid)
            v['midpoint_offsets_A'] = []
            v['midpoint_offsets_B'] = []
            for seg in v['A_segments']:
                v['midpoint_offsets_A'].append([0, 0, 0])
            for seg in v['B_segments']:
                v['midpoint_offsets_B'].append([0, 0, 0])
            
            
            if v['type'] =='PT':
                v['ghost_bisector'] = down
            elif v['type'] =='PB':
                v['ghost_bisector'] = up
            else: assert 0
            v['alpha'] = 0
            v['is_rigid'] = False
        
        if v['type'] == 'TH' or v['type'] == 'BH':
            # Ghost Normal already collected
            v['A_segments'] = []
            v['B_segments'] = []
            for eid in v['edges']:
                if edgeGlb[eid]['type'] == 'HT' or edgeGlb[eid]['type'] == 'NBH':
                    v['A_segments'].append(eid)
            v['midpoint_offsets_A'] = []
            v['midpoint_offsets_B'] = []
            for seg in v['A_segments']:
                v['midpoint_offsets_A'].append([0, 0, 0])
            if v['type'] =='TH':
                v['ghost_bisector'] = down # Not much reasoning here
            elif v['type'] =='BH':
                v['ghost_bisector'] = up # Not much reasoning here
            else: assert 0

            v['alpha'] = 0
            v['is_rigid'] = True
            

        if v['type'] == 'AA': # Arm joint
            v['A_segments'] = []
            v['B_segments'] = []
            ta_anb_pairs = []
            for eid_ta in v['edges']:
                if edgeGlb[eid_ta]['type'] == 'TA':
                    if vGlb[edgeGlb[eid_ta]['con'][0]]['type'] == 'PT':
                        v['ghost_normal'] = vGlb[edgeGlb[eid_ta]['con'][0]]['ghost_normal']
                    elif vGlb[edgeGlb[eid_ta]['con'][1]]['type'] == 'PT':
                        v['ghost_normal'] = vGlb[edgeGlb[eid_ta]['con'][1]]['ghost_normal']
                    else: assert 0
                    TA_vec = normalizeVec(get_edge_vec(vGlb, edgeGlb, eid_ta, vid))
                    v['aux_A'] = TA_vec
                    for eid_anb in v['edges']:
                        if edgeGlb[eid_anb]['type'] == 'ANB':
                            ANB_vec = normalizeVec(vecMult(get_edge_vec(vGlb, edgeGlb, eid_anb, vid), -1))
                            ta_anb_pairs.append([eid_ta, eid_anb])
            assert len(ta_anb_pairs) == 1
            v['B_segments'].append(ta_anb_pairs[0][0])
            v['B_segments'].append(ta_anb_pairs[0][1])

            v['midpoint_offsets_A'] = []
            v['midpoint_offsets_B'] = []
            for _ in v['B_segments']:
                v['midpoint_offsets_B'].append([0, 0, 0])
                    


            v['alpha'] = 0
            v['is_rigid'] = False
            v['ghost_bisector'] = down
            
            for eid in v['edges']:
                if edgeGlb[eid]['type'] == 'TA' or edgeGlb[eid]['type'] == 'ANB':
                    edgeGlb[eid]['segment_normal'] = normalizeVec(crossProduct3D(v['aux_A'], v['ghost_normal']))
        vGlb[vid] = v
    
    
    ## Define overhang segment normal as the same as the arm it is attached to
    for vid, v in enumerate(vGlb):
        if v['type'][0] == 'P':
            for eid_h in v['edges']:
                if edgeGlb[eid_h]['type'] == 'HT':
                    for eid_t in v['edges']:
                        if edgeGlb[eid_t]['type'] == 'TA':
                            edgeGlb[eid_h]['segment_normal'] = edgeGlb[eid_t]['segment_normal']
                            # print(edgeGlb[eid_h]['segment_normal'])
                if edgeGlb[eid_h]['type'] == 'NBH':
                    for eid_b in v['edges']:
                        if edgeGlb[eid_b]['type'] == 'ANB':
                            edgeGlb[eid_h]['segment_normal'] = edgeGlb[eid_b]['segment_normal']
                            # print(edgeGlb[eid_h]['segment_normal'])
                    
    return vGlb, edgeGlb

def handle_pivot_offsets(vGlb, edgeGlb, flipBit, thickness, overhang = None):
    # Always to be called after handle_offsets
    # Since it assures uid1 -> A segments, uid2 -> B segments
    # a FlipBit of 1 means, at rest state, the umbrella has the compliant joints closer to its central axis. 
    # That implies the AA joints must have an offset pointing inside. And the neighboring umbrellas have the offset pointing outside from them at the shared AA joint
    if overhang is not None: assert 0 # Not implemented

    for vid, v in enumerate(vGlb):
        if (v['type'] == 'PT' or v['type'] == 'PB') and (len(v['B_segments']) > 0):
            # Move inside towards plate
            
            seg_id = v['B_segments'][0]
            assert seg_id != -1 and len(v['B_segments']) == 1 and (edgeGlb[seg_id]['type'] == 'TA' or edgeGlb[seg_id]['type'] == 'ANB')
            
            ### TODO Possibly flip this chunk with the next TODO chunk. Investigate
            if v['type'] == 'PT':
                joint_disp = vecMult(normalizeVec(edgeGlb[seg_id]['segment_normal']), -thickness/2)
            else:
                joint_disp = vecMult(normalizeVec(edgeGlb[seg_id]['segment_normal']), thickness/2)
            ### TODO 
            
            vGlb[vid]['coord'] = vecAdd(vGlb[vid]['coord'], joint_disp)
            pivot_offset = vecMult(joint_disp, -1)
            for seg_id in range(len(vGlb[vid]['midpoint_offsets_B'])):
                vGlb[vid]['midpoint_offsets_B'][seg_id] = vecAdd(vGlb[vid]['midpoint_offsets_B'][seg_id], pivot_offset)
            # Lower or lift depending on top or bottom plate
            seg_id = v['A_segments'][0]
            assert len(v['A_segments']) == 1 and (edgeGlb[seg_id]['type'] == 'PT' or edgeGlb[seg_id]['type'] == 'PB')
            
            joint_disp = vecMult(normalizeVec(edgeGlb[seg_id]['segment_normal']), -thickness/2)
            
            ### TODO 
            vGlb[vid]['coord'] = vecAdd(vGlb[vid]['coord'], joint_disp)
            ### TODO
            
            pivot_offset = vecMult(joint_disp, -1)
            for seg_id in range(len(vGlb[vid]['midpoint_offsets_A'])):
                vGlb[vid]['midpoint_offsets_A'][seg_id] = vecAdd(vGlb[vid]['midpoint_offsets_A'][seg_id], pivot_offset)
            
        if v['type'] == 'AAM':
            uid1, uid2 = v['uid']
            if flipBit[uid1]:
                seg_ids = v['A_segments']
            else:
                assert flipBit[uid2]
                seg_ids = v['B_segments']
            assert edgeGlb[seg_ids[0]]['segment_normal'] == edgeGlb[seg_ids[1]]['segment_normal']
            
            joint_disp = vecMult(edgeGlb[seg_ids[0]]['segment_normal'], -thickness/2)    
            vGlb[vid]['coord'] = vecAdd(vGlb[vid]['coord'], joint_disp)
            pivot_offset = vecMult(joint_disp, -1)        
            for seg_id in range(len(vGlb[vid]['midpoint_offsets_A'])):
                vGlb[vid]['midpoint_offsets_A'][seg_id] = vecAdd(vGlb[vid]['midpoint_offsets_A'][seg_id], pivot_offset)
            for seg_id in range(len(vGlb[vid]['midpoint_offsets_B'])):
                vGlb[vid]['midpoint_offsets_B'][seg_id] = vecAdd(vGlb[vid]['midpoint_offsets_B'][seg_id], pivot_offset)
    return vGlb
    
def handle_offsets(vGlb, edgeGlb, nbrMap, armJointAxisOffset):
    for joint1_id in range(len(vGlb)):
        if (vGlb[joint1_id]['type'] != 'AA'):
            continue
        uid1 = vGlb[joint1_id]['uid'][0]
        for joint2_id in range(len(vGlb)):
            
            if joint1_id == joint2_id: continue
            if (vGlb[joint2_id]['type'] != 'AA'):
                continue
            uid2 = vGlb[joint2_id]['uid'][0]
            if uid2 not in nbrMap[uid1]:
                continue
            if uid2 != vGlb[joint1_id]['uid'][1] or uid1 != vGlb[joint2_id]['uid'][1]: continue
            if abs(vecLength(vecSub(vGlb[joint1_id]['coord'], vGlb[joint2_id]['coord'])) - 2*armJointAxisOffset) > 1e-8: 
                continue
            
            assert uid2 == vGlb[joint1_id]['uid'][1]
            assert uid1 == vGlb[joint2_id]['uid'][1]

            
            
            #Merge these two joints
            vGlb[joint1_id]['type'] = 'AAM' # AA Merged
            vGlb[joint1_id]['coord'] = vecMult(vecAdd(vGlb[joint1_id]['coord'], vGlb[joint2_id]['coord']), 0.5)
            for eid in vGlb[joint2_id]['B_segments']:
                if edgeGlb[eid]['con'][0] == joint2_id:
                    edgeGlb[eid]['con'][0] = joint1_id
                elif edgeGlb[eid]['con'][1] == joint2_id:
                    edgeGlb[eid]['con'][1] = joint1_id
                
                    
            vGlb[joint1_id]['A_segments'] = copy.deepcopy(vGlb[joint2_id]['B_segments'])
            
            # A segments correspond to first uid
            vGlb[joint1_id]['uid'].reverse()
            assert vGlb[joint1_id]['uid'][0] == edgeGlb[vGlb[joint1_id]['A_segments'][0]]['uid']
            assert vGlb[joint1_id]['uid'][1] == edgeGlb[vGlb[joint1_id]['B_segments'][0]]['uid']

            vGlb[joint1_id]['midpoint_offsets_A'] = []
            vGlb[joint1_id]['midpoint_offsets_B'] = []
            for _ in vGlb[joint1_id]['A_segments']:
                vGlb[joint1_id]['midpoint_offsets_A'].append(vecSub(vGlb[joint2_id]['coord'], vGlb[joint1_id]['coord']))
            for _ in vGlb[joint1_id]['B_segments']:
                vGlb[joint1_id]['midpoint_offsets_B'].append(vecMult(vGlb[joint1_id]['midpoint_offsets_A'][0], -1.0))
    newVGlb = []
    new_vid = [-1]*len(vGlb)
    for old_vid in range(len(vGlb)):
        if vGlb[old_vid]['type'] == 'AA': # remove
            pass
        else:
            new_vid[old_vid] = len(newVGlb)
            newVGlb.append(copy.deepcopy(vGlb[old_vid]))
            
    
    newEdgeGlb = []
    new_eid = [-1]*len(edgeGlb)
    for old_eid in range(len(edgeGlb)):
        ov1, ov2 = edgeGlb[old_eid]['con']
        edgeGlb[old_eid]['con'] = [new_vid[ov1], new_vid[ov2]]
        if new_vid[ov1] == -1 or new_vid[ov2] == -1: # Remove
            assert 0
            pass
        else:
            new_eid[old_eid] = len(newEdgeGlb)
            newEdgeGlb.append(edgeGlb[old_eid])
    for vid in range(len(newVGlb)):
        for eid in range(len(newVGlb[vid]['edges'])):
            newVGlb[vid]['edges'][eid] = new_eid[newVGlb[vid]['edges'][eid]]
        for seg_id in range(len(newVGlb[vid]['A_segments'])):
            newVGlb[vid]['A_segments'][seg_id] = new_eid[newVGlb[vid]['A_segments'][seg_id]]
        for seg_id in range(len(newVGlb[vid]['B_segments'])):
            newVGlb[vid]['B_segments'][seg_id] = new_eid[newVGlb[vid]['B_segments'][seg_id]]

    for v in newVGlb:
        assert len(v['A_segments']) == len(v['midpoint_offsets_A'])
        assert len(v['B_segments']) == len(v['midpoint_offsets_B']), v['type']
    
    return newVGlb, newEdgeGlb

    
def get_bbox_diagonal(vertices):
    bmin = copy.copy(vertices[0])
    bmax = copy.copy(vertices[0])
    for vert in vertices:
        for i in range(3):
            bmin[i] = min(vert[i], bmin[i])
            bmax[i] = max(vert[i], bmax[i])
    return vecLength(vecSub(bmax, bmin))

def scale_vecs(vecs, scale):
    if isinstance(vecs, list):
        if None in vecs:
            for vi, vec in enumerate(vecs):
                if vec is not None:
                    vecs[vi] = vecMult(vec, scale)
        else:
            if len(vecs) == 0:
                return vecs
            if isinstance(vecs[0], list):
                for vi, vec in enumerate(vecs):
                    if len(vec) == 0:
                        continue
                    if isinstance(vec[0], list):
                        for vj, vec in enumerate(vecs[vi]):
                            vecs[vi][vj] = vecMult(vecs[vi][vj], scale)
                    else:
                        vecs[vi] = vecMult(vec, scale)
            else:
                if None in vecs:
                    for vi, vec in enumerate(vecs):
                        if vec is not None:
                            vecs[vi] = vecMult(vec, scale)
                else:
                    vecs = vecMult(vecs, scale)
    else:
        vecs = scale * vecs
    return vecs

def scale_data(input_data, scale):
    input_data['vertices'] = scale_vecs(input_data['vertices'], scale)
    input_data['midpoint_offsets_A'] = scale_vecs(input_data['midpoint_offsets_A'], scale)
    input_data['midpoint_offsets_B'] = scale_vecs(input_data['midpoint_offsets_B'], scale)
    input_data['correspondence'] = scale_vecs(input_data['correspondence'], scale)
    input_data['plate_correspondence'] = scale_vecs(input_data['plate_correspondence'], scale)
    input_data['boundary_correspondence'] = scale_vecs(input_data['boundary_correspondence'], scale)
    input_data['base_mesh_v'] = scale_vecs(input_data['base_mesh_v'], scale)
    input_data['target_v'] = scale_vecs(input_data['target_v'], scale)
    input_data['plate_edge_length'] = scale_vecs(input_data['plate_edge_length'], scale)
    input_data['arm_plate_edge_offset'] = scale_vecs(input_data['arm_plate_edge_offset'], scale)
    input_data['arm_joint_offset'] = scale_vecs(input_data['arm_joint_offset'], scale)
    input_data['margin_length'] = scale_vecs(input_data['margin_length'], scale)
    input_data['thickness'] = scale_vecs(input_data['thickness'], scale)
    input_data['width'] = scale_vecs(input_data['width'], scale)
    return input_data

def export_data(vGlb, edgeGlb, flipBit, t_mesh, uv_mesh, scaleFactors, nbr_map, v_out, f_out, plateLength, armPlateEdgeAxisOffset, armJointAxisOffset, width, thickness, targetSpacingFactor, marginLength, filename = 'data.json', select_umbrella = None, min_coerce_dist = 1e-8):
    input_data = {'vertices': [],
                  'edges' : [],
                  'alphas' : [],
                  'ghost_bisectors' : [],
                  'ghost_normals' : [],
                  'A_segments' : [],
                  'B_segments' : [],
                  'midpoint_offsets_A' : [],
                  'midpoint_offsets_B' : [],
                  'segment_normals' : [],
                  'is_rigid': [],
                  'uid': [],
                  'v_labels': [],
                  'e_labels': [],
                  'scale_factors': [],
                  'correspondence': [],
                  'plate_correspondence' : [[] for _ in range(len(nbr_map))],
                  'boundary_correspondence' : [None for _ in range(len(nbr_map))],
                  'thickness': thickness,
                  'width': width,
                  'target_spacing_factor': targetSpacingFactor}
                            
    if select_umbrella == True:
        vGlb_new = []
        new_vid = [-1 for _ in range(len(vGlb))]
        for old_vid, v in enumerate(vGlb):
            if isinstance(v['uid'], list):
                if select_umbrella in v['uid']:
                    vGlb_new.append(copy.deepcopy(v))
                    new_vid[old_vid] = len(vGlb_new) - 1
            if v['uid'] == select_umbrella:
                vGlb_new.append(copy.deepcopy(v))
                new_vid[old_vid] = len(vGlb_new) - 1
        edgeGlb_new = []
        new_eid = [-1 for _ in range(len(edgeGlb))]
        for old_eid, e in enumerate(edgeGlb):
            vid1, vid2 = e['con']
            if new_vid[vid1] == -1 or new_vid[vid2] == -1:
                pass
            else:
                e['con'] = [new_vid[vid1], new_vid[vid2]]
                edgeGlb_new.append(copy.deepcopy(e))
                new_eid[old_eid] = len(edgeGlb_new) - 1
        for vid in range(len(vGlb_new)):
            for seg_id in range(len(vGlb_new[vid]['A_segments'])):
                vGlb_new[vid]['A_segments'][seg_id] = new_eid[vGlb_new[vid]['A_segments'][seg_id]]
            for seg_id in range(len(vGlb_new[vid]['B_segments'])):
                vGlb_new[vid]['B_segments'][seg_id] = new_eid[vGlb_new[vid]['B_segments'][seg_id]]
        for v in vGlb_new:
            for seg_id in range(len(v['A_segments'])-1, -1, -1):
                if v['A_segments'][seg_id] == -1:
                    v['A_segments'].remove(v['A_segments'][seg_id])
                    v['midpoint_offsets_A'].remove(v['midpoint_offsets_A'][seg_id])
            for seg_id in range(len(v['B_segments'])):
                if v['B_segments'][seg_id] == -1:
                    v['B_segments'].remove(v['B_segments'][seg_id])
                    v['midpoint_offsets_B'].remove(v['midpoint_offsets_B'][seg_id])
            v['uid'] = 0
        v_out_new = [] 
        for vid, v in enumerate(v_out):
            if vid in f_out[select_umbrella]:
                v_out_new.append(v)
        assert len(v_out_new) == 3
    
        f_out = [[0,1,2]]
        nbr_map = [[]]
        vGlb = vGlb_new
        edgeGlb = edgeGlb_new
        v_out = v_out_new



    # Global config Data
    input_data['plate_edge_length'] = plateLength
    input_data['arm_plate_edge_offset'] = armPlateEdgeAxisOffset
    input_data['arm_joint_offset'] = armJointAxisOffset
    input_data['margin_length'] = marginLength
    input_data['scale_factors'] = scaleFactors

    # Sim Data
    #   
    for vid, _ in enumerate(vGlb):
        v = copy.deepcopy(vGlb[vid])
        input_data['vertices'].append(v['coord'])
        input_data['alphas'].append(v['alpha'])
        input_data['ghost_bisectors'].append(v['ghost_bisector'])
        input_data['ghost_normals'].append(v['ghost_normal'])
        input_data['A_segments'].append(v['A_segments'])
        input_data['B_segments'].append(v['B_segments'])
        input_data['midpoint_offsets_A'].append(v['midpoint_offsets_A'])
        input_data['midpoint_offsets_B'].append(v['midpoint_offsets_B'])
        input_data['is_rigid'].append(v['is_rigid'])
        input_data['uid'].append(v['uid'])
        input_data['correspondence'].append(None)

    for e in edgeGlb:
        input_data['edges'].append(e['con'])
        input_data['segment_normals'].append(e['segment_normal'])

    # Vis Data
    for v in vGlb:
        input_data['v_labels'].append(v['type'])
    for e in edgeGlb:
        input_data['e_labels'].append(e['type'])
    input_data['flip_bits'] = flipBit
    input_data['base_mesh_v'] = v_out
    input_data['base_mesh_f'] = f_out
    input_data['target_v'] = []
    input_data['target_f'] = []

    if t_mesh is not None:
        target_v = t_mesh.Vertices.ToPoint3dArray()
        target_f = t_mesh.Faces
        for t_v in target_v:
            input_data['target_v'].append([t_v[0], t_v[1], t_v[2]])
        for t_fid in range(target_f.Count):
            t_f = target_f.Item[t_fid]
            input_data['target_f'].append([t_f[0], t_f[1], t_f[2]])

    input_data['bbox_diagonal'] = get_bbox_diagonal(input_data['vertices'])
    scale_down = 1.0/input_data['bbox_diagonal']
    input_data = scale_data(input_data, scale_down)

    
    with gzip.open(filename, 'wt') as outfile:
        json.dump(input_data, outfile, indent=4)

    vGlb_scaled = copy.deepcopy(vGlb)
    for vid, v in enumerate(vGlb_scaled):
        vGlb_scaled[vid]['coord'] = scale_vecs(vGlb_scaled[vid]['coord'], scale_down)
        vGlb_scaled[vid]['midpoint_offsets_A'] = scale_vecs(vGlb_scaled[vid]['midpoint_offsets_A'], scale_down)
        vGlb_scaled[vid]['midpoint_offsets_B'] = scale_vecs(vGlb_scaled[vid]['midpoint_offsets_B'], scale_down)

    return vGlb, edgeGlb, vGlb_scaled
def get_checkerboard(nbr_map, seed_uid = 0):
    num_uids = len(nbr_map)
    # Since the pivot is asymmetric, the cells form a checkerboard pattern
    flipBit = [False for uid in range(num_uids)]
    visitBit = [False for uid in range(num_uids)]
    flipBit[seed_uid] = True
    visitBit[seed_uid] = True
    while False in visitBit:
        for uid1 in range(num_uids):
            if not visitBit[uid1]: continue
            for uid2 in nbr_map[uid1]:
                if flipBit[uid2] == flipBit[uid1]:
                    flipBit[uid2] = not flipBit[uid1]
                visitBit[uid2] = True
    return flipBit


                
def heightFromConformalScale(scale, thickness, targetSpacingFactor, plateLength, marginLength, flip, handlePivots = True, degree = 3):
    assert scale > 1.0

    if degree == 3:
        height = (scale - 1.0)*plateLength/sqrt(3.) + 2*marginLength*scale
    elif degree == 4:
        height = (scale - 1.0)*plateLength + 2*marginLength*scale
    elif degree == 6:
        height = (scale - 1.0)*plateLength*sqrt(3.) + 2*marginLength*scale
    
    # if not flip: #thickness footprint is accounted in this umbrella
    #     height +=  2 * scale * thickness
    
    # If splitting accountability for thickness symmetrically despite pivot joint asymmetry
    if handlePivots:
        height +=  scale * thickness

    height = heightUnProj(height, targetSpacingFactor, thickness)
    # pivot offset handling would move the TA joints of neighboring plates inward so that a spacing of "thickness" is made between them. Similarly for the z direction
    if handlePivots:
        height+= thickness
    return height

def getEdgeLength(plateLength, marginLength, thickness, handlePivots = True, degree = 3):
    if degree == 3:
        edgeLength = plateLength + 2*sqrt(3)*marginLength
        if handlePivots: edgeLength += 2*sqrt(3)*thickness/2
    elif degree == 4:
        edgeLength = plateLength + 2*marginLength
        if handlePivots: edgeLength += thickness
    elif degree == 6:
        edgeLength = plateLength + 2*marginLength/sqrt(3)
        if handlePivots: edgeLength += thickness/sqrt(3)
    return edgeLength


def getPlateLength(edgeLength,
                   marginLength,
                   thickness,
                   handlePivots = True,
                   degree = 3):
    if degree == 3:
        plateLength = edgeLength - 2*sqrt(3.)*marginLength # [RQ]
        if handlePivots: plateLength -= 2*sqrt(3)*thickness/2
    elif degree == 4:
        plateLength = edgeLength - 2*marginLength
        if handlePivots: plateLength -= thickness
    elif degree == 6:
        plateLength = edgeLength - 2*marginLength/sqrt(3)
        if handlePivots: plateLength -= thickness/sqrt(3)
    return plateLength

def heightProj(height, targetSpacingFactor, thickness):
    s = targetSpacingFactor*thickness
    return sqrt(height**2 - (s - thickness)**2)

def heightUnProj(height, targetSpacingFactor, thickness):
    s = targetSpacingFactor*thickness
    return sqrt(height**2 + (s - thickness)**2)

def minHeightWithoutScaling(scaleFactors, plateLength, thickness, targetSpacingFactor, handlePivots = True, degree = 3):
    minScaleFactor = min(scaleFactors)
    if minScaleFactor <= 1: # Without uniform scaling, scaleFactors should be greater than 1 as umbrellas can only expand
            return -1
    if degree == 3:
        if not handlePivots:
            return plateLength*(minScaleFactor - 1.0)/sqrt(3.)
        return  plateLength*(minScaleFactor - 1.0)/sqrt(3.) + minScaleFactor * thickness  
    elif degree == 4:
        if not handlePivots:
            return plateLength*(minScaleFactor - 1.0)
        return  plateLength*(minScaleFactor - 1.0) + minScaleFactor * thickness  
    elif degree == 6:
        if not handlePivots:
            return plateLength*(minScaleFactor - 1.0)*sqrt(3.)
        return  plateLength*(minScaleFactor - 1.0)*sqrt(3.) + minScaleFactor * thickness  

def getTargetScaling(scaleFactors, minHeight, plateLength, thickness, targetSpacingFactor, handlePivots = True, degree = 3):
    minScaleFactor = min(scaleFactors)
    if degree == 3:
        if not handlePivots:
            return (plateLength + sqrt(3)*heightProj(minHeight, targetSpacingFactor, thickness))/(plateLength)/minScaleFactor
        # If splitting accountability for thickness symmetrically despite pivot joint asymmetry
        return (plateLength + sqrt(3)*heightProj(minHeight, targetSpacingFactor, thickness))/(plateLength + sqrt(3)*thickness)/minScaleFactor
        # # If splitting accountability for thickness asymmetrically because of pivot joint asymmetry
        # return (1 + sqrt(3)*minHeight/plateLength)/minScaleFactor
    elif degree == 4:
        if not handlePivots:
            return (plateLength + heightProj(minHeight, targetSpacingFactor, thickness))/(plateLength)/minScaleFactor
        # If splitting accountability for thickness symmetrically despite pivot joint asymmetry
        return (plateLength + heightProj(minHeight, targetSpacingFactor, thickness))/(plateLength + thickness)/minScaleFactor
        # # If splitting accountability for thickness asymmetrically because of pivot joint asymmetry
        # return (1 + minHeight/plateLength)/minScaleFactor
    elif degree == 6:
        if not handlePivots:
            return (plateLength + (1/sqrt(3))*heightProj(minHeight, targetSpacingFactor, thickness))/(plateLength)/minScaleFactor
        # If splitting accountability for thickness symmetrically despite pivot joint asymmetry
        return (plateLength + (1/sqrt(3))*heightProj(minHeight, targetSpacingFactor, thickness))/(plateLength + thickness/sqrt(3))/minScaleFactor
        # # If splitting accountability for thickness asymmetrically because of pivot joint asymmetry
        # return (1 + (1/sqrt(3))*minHeight/plateLength)/minScaleFactor

def getVertexNbrMap(faces, num_vertices):
    nbr_map = [[] for i in range(num_vertices)]
    for face in faces:
        for vi in face:
            for vj in face:
                if vi != vj and vj not in nbr_map[vi]:
                    nbr_map[vi].append(vj)
    return nbr_map
def getScaleFactorPerVertex(scale_factors, faces, num_vertices):
    per_vertex = [0 for i in range(num_vertices)]
    counts = [0 for i in range(num_vertices)]
    for fid, face in enumerate(faces):
        for vi in face:
            per_vertex[vi] += scale_factors[fid]
            counts[vi] += 1
    for vi in range(len(per_vertex)):
        per_vertex[vi] /= counts[vi]
    return per_vertex




def genPattern(edgeLength,
               t_mesh,
               uv_mesh,
               i_out,
               v_out,
               f_out,
               c_out,
               x_out,
               scaleLength,
               jsonPath,
               marginLength = 0.0,
               armPlateEdgeAxisOffset = 0.0,
               armJointAxisOffset = 0.0,
               asymmetryOffset = 0,
               width = 5,
               thickness = 3,
               targetSpacingFactor = 1,
               minHeight = 80,
               select_umbrella = False,
               handlePivots = True,
               min_coerce_dist = 1e-8,
               degree = 3,
               overhang = None,
               heights = None):

    mesh_degree = degree
    # pivot offset handling would move the TA joints of neighboring plates inward so that a spacing of "thickness" is made between them. Similarly for the z direction
    plateLength = getPlateLength(edgeLength, marginLength, thickness, handlePivots, degree=mesh_degree)
    
    vList = []
    topList = []
    v_count_gt = 0
    
    nbr_map = copy.copy(x_out)
    
    # print("min_id: " + str(scaleLength.index(min(scaleLength))))
    # MinHeight umbrella must have flipBit = 1, => thickness footprint is not accounted in this umbrella => height of this umbrella wouldn't get a 2*s*t factor added in heightFromConformalScale
    flipBit = get_checkerboard(nbr_map, seed_uid=scaleLength.index(min(scaleLength)))
    if minHeight == 0: # no uniform conformal scaling
        pass
    else:
        # # If splitting accountability for thickness asymmetrically because of pivot joint asymmetry
        # scaleLength = [x*getTargetScaling(scaleLength, minHeight, plateLength, thickness, handlePivots) for x in scaleLength]

        # If splitting accountability for thickness symmetrically despite pivot joint asymmetry
        scaleLength = [x*getTargetScaling(scaleLength, minHeight, plateLength, thickness, targetSpacingFactor, handlePivots) for x in scaleLength]
    
    # Transform to dual for hex topology
    v_nbr_map = getVertexNbrMap(f_out, len(v_out))
        
    if degree == 3 or degree == 4:
        data = c_out
        scale_data = scaleLength
    elif degree == 6:
        data = v_out
        if len(scaleLength) == len(f_out):
            scale_data = getScaleFactorPerVertex(scaleLength, f_out, len(v_out))
        elif len(scaleLength) == len(v_out):
            scale_data = scaleLength
        else: assert 0
        nbr_map = v_nbr_map
        flipBit = get_checkerboard(v_nbr_map, seed_uid=scale_data.index(min(scale_data)))
    for j,(c, scale_factor) in enumerate(zip(data, scale_data)):
        v_count_gt += 2 + mesh_degree * 2
        v_count_gt += len(nbr_map[j]) * (1 + 2 * (overhang is not None))
        
        if heights is None:
            zHeight = heightFromConformalScale(scale_factor, thickness, targetSpacingFactor, plateLength, marginLength, flip = flipBit[j], handlePivots=handlePivots)
        else:
            zHeight = heights[j]
        zHeight_nb = []
        c_nb = []
        # print(zHeight, overhang[j], zHeight + overhang[j])
        for nb_id in range(len(nbr_map[j])):
            if heights is None:
                zHeight_nb.append(heightFromConformalScale(scale_data[nbr_map[j][nb_id]], thickness, targetSpacingFactor, plateLength, marginLength, flip = flipBit[nbr_map[j][nb_id]], handlePivots=handlePivots))
            else:
                zHeight_nb.append(heights[nbr_map[j][nb_id]])
            c_nb.append(data[nbr_map[j][nb_id]])
        # generate all the umbrellae
        
        vL, topL = genLinkageVerts(zHeight, zHeight_nb, j, nbr_map, c, c_nb, edgeLength, marginLength, armPlateEdgeAxisOffset, armJointAxisOffset, asymmetryOffset, degree = mesh_degree, overhang=overhang)
        
        vList.append(vL)
        topList.append(topL)
        
            
    [rsEdge, edgeList] = genLinkage(vList, topList, degree = mesh_degree, overhang=overhang)
    
    vGlb, edgeGlb, _ = collect_all(vList, edgeList)
    assert len(vGlb) == v_count_gt
    
    vGlb, edgeGlb = collect_ghost_data(vGlb, edgeGlb, degree = mesh_degree)
    vGlb, edgeGlb = handle_offsets(vGlb, edgeGlb, nbr_map, armJointAxisOffset)
    if handlePivots:
        vGlb = handle_pivot_offsets(vGlb, edgeGlb, flipBit, thickness, overhang=overhang)

    vGlb_exp, edgeGlb_exp, vGlb_scaled = export_data(vGlb, edgeGlb, flipBit, t_mesh, uv_mesh, scaleLength, nbr_map, v_out, f_out, plateLength, armPlateEdgeAxisOffset, armJointAxisOffset, width, thickness, targetSpacingFactor, marginLength, filename=jsonPath, select_umbrella=select_umbrella, min_coerce_dist=min_coerce_dist)
    assert vGlb == vGlb_exp and edgeGlb == edgeGlb_exp
    
    
    return [vGlb_exp, edgeGlb]


