from math import sqrt, tan, sin, cos, atan, asin, acos, pi, atan2
def flatten(l):
    return [item for sublist in l for item in sublist]
    
def vecAdd(v1,v2):
    return [x+y for x,y in zip(v1,v2)]

def vecSub(v1,v2):
    return [x-y for x,y in zip(v1,v2)]

def vecMult(v1,a):
    return [x*a for x in v1]
    
def vecLength (v):
    return sqrt(sum([x**2 for x in v]))

def dotProduct(v1, v2):
    prod = [x*y for x,y in zip(v1,v2)]
    return sum([x for x in prod])
def crossProduct3D(v1, v2):
    assert len(v1) == 3 and len(v2) == 3
    a1, a2, a3 = v1
    b1, b2, b3 = v2
    return [a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1]

def angleBWVecs(v1, v2):
    return acos(dotProduct(v1, v2)/(vecLength(v1)*vecLength(v2)))

def normalizeVec(v1):
    return [x/vecLength(v1) for x in v1]

def vecBisector(v1, v2):
    if round(vecLength(v1), 8) != 1 or round(vecLength(v2), 8) != 1:
        print("Warning!: vecBisector called with unnormalized vecs. Normalizing")
    v1_norm = normalizeVec(v1)
    v2_norm = normalizeVec(v2)
    return normalizeVec(vecAdd(v1_norm, v2_norm))

def vDic(v):
    return {'type': v[0], 'edges': v[1], 'coord': v[2]}
    
def ptRot(v,a):
    return vDic([v['type'], v['edges'],  [v['coord'][0]*cos(a) - v['coord'][1]*sin(a), v['coord'][0]*sin(a) + v['coord'][1]*cos(a), v['coord'][2]]] )
    
def ptRotL(vL,a):
    return [ptRot(v,a) for v in vL]
    
def ptRotTri(v):
    vL = []
    vL.append(v)
    vL.append(ptRot(v,2.*pi/3))
    vL.append(ptRot(v,4.*pi/3))
    return vL

def ptRotPoly(v, degree = 3):
    vL = []
    vL.append(v)
    for i in range(degree-1):
        vL.append(ptRot(v, 2.0 * (i+1) * pi/degree))
    return vL
    
def ppClosest(h, n, edgeLength):
    """
    Find the point in h that n is closest to
    """
    tempD = [vvDist(hi,n) for hi in h]
    if min(tempD)>edgeLength/10.:
        errorTog = True
    else:
        errorTog = False
    return tempD.index(min(tempD)), errorTog

# creates a lever in the XZ plane that is moved to start and rotated
def vvDist(v1,v2):
    return sqrt(sum([(v2['coord'][i]-v1['coord'][i])**2 for i in [0,1,2]]))
   
# creates a lever in the XZ plane that is moved to start and rotated
def ppDist(v1,v2):
    return sqrt(sum([(v2[i]-v1[i])**2 for i in [0,1,2]]))
   
def list_to_dict(a):
    it = iter(a)
    res_dct = dict(zip(it, it))
    return res_dct


def moveV(v,a):
    return vDic([v['type'], v['edges'], [v['coord'][i]+a[i] for i in [0,1,2]]])

def mirrV(v,c,plane):
    if plane == 'XZ':
        return  vDic([v['type'], v['edges'], [v['coord'][0], v['coord'][1]+2*(c[1]-v['coord'][1]), v['coord'][2]]])
    else:
        assert 0

def avgPt(p):
    """
    Find the centroid of a given list of points
    """
    x=[ sum([p_i[i] for p_i in p])/len(p) for i in [0,1,2]]
    
    return x
    
def noneList(x,a):
    """ change the indices a of x to [-1] """
    for ai in a:
        x[ai] = [-1]
        
def noneListR(a):
    """ remove [-1] """
    idxD = dict()
    idx = 0
    temp = []
    for i, ai in enumerate(a):
        if ai != [-1]:
            temp.append(ai)
            idxD[i]=idx
            idx+=1
            
    return temp, idxD

def generate_regular_sqr_grid(rows, cols, l=1, ox=0, oy=0):
    
    v_grid = [] # vertex list of x,y,z coord
    f_grid = [] # per faces has lists of vertice ids
    c_grid = [] # centroid coordinate per element
    i_grid = [] # due to the mirroring nature of the aux and bass unit cels, we record whether 1) tri, up or down, 2) quad, top left, bottom left, top right, bottom right
    x_grid = [] # per face, the id of all the neighbouring faces
    
    h = l # Height of the square normalized to 1
    
    # generate all the vertices for the squares
    idE = 0 # Counter for elements
    for row in range(rows):
        for col in range(cols):
            v_grid.append([l*col+ox,h*row+oy, 0])
    
    # generate all the connectivities (node numbers)
    # and the quadrant per square
    for row in range(rows-1):
        for col in range(cols-1):
            
            n = [col + row*cols,
                 col + row*cols + 1,
                 col + (row+1)*cols + 1,
                 col + (row+1)*cols]
            
            if row%2 == 0:
                if col%2 == 0:
                    fL = [n[0], n[1], n[2], n[3]]
                    i_grid.append([ 1, 1, -1,  -1])
                else:
                    fL = [n[1], n[0], n[3], n[2]]
                    i_grid.append([-1, 1,  1,  -1])
            else:
                if col%2 == 0:
                    fL = [n[3], n[2], n[1], n[0]]
                    i_grid.append([ 1, -1, -1,  1])
                else:
                    fL = [n[2], n[3], n[0], n[1]]
                    i_grid.append([-1, -1,  1,  1])
            
            # cL is the centroid per square
            cL = avgPt([v_grid[fLi] for fLi in fL])
            
            f_grid.append(fL)
            c_grid.append(cL)
                
            # Find the ID of the neighbouring elements (special cases are the edges, and corners)
            if row == 0:
                if col == 0:
                    x_grid.append([                          idE + 1, idE + cols - 1])
                elif col == cols-2:
                    x_grid.append([idE - 1,                           idE + cols - 1])
                else:
                    x_grid.append([idE - 1,                  idE + 1, idE + cols - 1])
            elif row == rows-2:
                if col == 0:
                    x_grid.append([          idE - cols + 1, idE + 1                ])
                elif col == cols-2:
                    x_grid.append([ idE - 1, idE - cols + 1                         ])
                else:
                    x_grid.append([ idE - 1, idE - cols + 1, idE + 1                ])
            else:
                if col == 0:
                    x_grid.append([          idE - cols + 1, idE + 1, idE + cols - 1])
                elif col == cols-2:
                    x_grid.append([ idE - 1, idE - cols + 1,          idE + cols - 1])
                else:
                    x_grid.append([ idE - 1, idE - cols + 1, idE + 1, idE + cols - 1])

            idE+=1   

    return v_grid, f_grid, c_grid, i_grid, x_grid

def generate_regular_tri_grid(rows, cols, l=1, ox=0, oy=0):
    
    v_grid = [] # vertex list of x,y,z coord
    f_grid = [] # per faces has lists of vertice ids
    c_grid = [] # centroid coordinate per element
    i_grid = [] # due to the mirroring nature of the aux and bass unit cels, we record whether 1) tri, up or down, 2) quad, top left, bottom left, top right, bottom right
    x_grid = [] # per face, the id of all the neighbouring faces
    
    h = sqrt(3)/2*l # Height of the triangle
    
    idE = 0 # Counter for elements
    for row in range(rows):
        for col in range(cols):
            if row%2==0:
                v_grid.append([l*col+ox,h*row+oy, 0])
            else:
                v_grid.append([l*col+l/2+ox,h*row+oy, 0])
                
    for row in range(rows-1):
        iList = []
        for col in range(cols-1):
            
            if row%2==0:
                m = 1
                n = 0
                iList.append(0)
                iList.append(1)
                
            else:
                m = 0
                n = 1
                iList.append(1)
                iList.append(0)
        
            n0 = col + row*cols
            n1 = col + row*cols + 1
            n2 = col + (row+1)*cols + n
            
            # iso triangle #1
            fL1 = [n0,n1,n2]
        
            c0 = (v_grid[n0][0]+v_grid[n1][0])/2
            c1 = v_grid[n0][1] + sqrt(3)/6*l
            
            # center of tri #1
            cL1 = [c0,c1,0]
            
            n0 = col + row*cols + m
            n1 = col + (row+1)*cols + 1
            n2 = col + (row+1)*cols
            
            # iso triangle #2
            fL2 = [n0,n1,n2]
            
            c0 = (v_grid[n1][0]+v_grid[n2][0])/2
            c1 = v_grid[n1][1] - sqrt(3)/6*l
            
            # center of tri #2
            cL2 = [c0,c1,0]

            # this generation always starts with an a triangle with bottom edge aligned with the xaxis first
            if row%2==0:
                    
                f_grid.append(fL1)
                c_grid.append(cL1)
                f_grid.append(fL2)
                c_grid.append(cL2)
                
                # Find the ID of the neighbouring elements
                if row == 0:
                    if col == 0:
                        x_grid.append([idE+1])
                    else:
                        x_grid.append([idE-1, idE+1])
                else:
                    if col ==0:
                        x_grid.append([idE+1, idE-(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE-(cols-1)*2, idE+1])
                     
                idE+=1   
                
                # Find the ID of the neighbouring elements
                if row == rows-2:
                    if col == cols-2:
                        x_grid.append([idE-1])
                    else:
                        x_grid.append([idE-1, idE+1])
                else:
                    if col == cols-2:
                        x_grid.append([idE-1, idE+(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE+(cols-1)*2, idE+1])
                        
                idE+=1
                
            else:
                f_grid.append(fL2)
                c_grid.append(cL2)
                
                f_grid.append(fL1)
                c_grid.append(cL1)
                
                # Find the ID of the neighbouring elements
                if row == rows - 2:
                    if col ==0:
                        x_grid.append([idE+1])
                    else:
                        x_grid.append([idE-1, idE+1])

                else:
                    if col ==0:
                        x_grid.append([idE+1, idE+(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE+(cols-1)*2, idE+1])
                    
                idE+=1
                
                # Find the ID of the neighbouring elements
                if row == 0:
                    # But this wont happen unless we flip the order or even and odd
                    if col == cols - 2:
                        x_grid.append([idE-1])
                    else:
                        x_grid.append([idE-1, idE+1])
                else:
                    if col == cols-2:
                        x_grid.append([idE-1, idE-(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE-(cols-1)*2, idE+1])
                        
                idE+=1
                
        # in even last one, i dont exclude the upper triangle, which causes issues
        i_grid.append(iList)
        
    i_grid = [item for sublist in i_grid for item in sublist]
    
    return v_grid,f_grid,c_grid,i_grid, x_grid