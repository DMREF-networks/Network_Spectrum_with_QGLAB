import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import random

####Notes:
###latmat is a dxd matrix whose columns are lattice vectors

###################################################################################################
### Functions to make point patterns (2D)

## URL (Z2)
#Stay less than or roughly a = 0.3 to make sure the edges are not too short
def gen_URL(N, a):
    sites = np.floor(np.sqrt(N)).astype(int) 
    coords = np.empty([sites**2,2])
    for i in range(coords.shape[0]):
        coords[i] = [i//sites,i%sites]

    if a != 0:
        for i in range(coords.shape[0]):
            coords[i] += (np.random.rand(2) - 0.5) * a

    latmat = np.array([[sites,0],[0,sites]])
    wrap(coords, latmat)

    return(coords, latmat)

def gen_A2(N,a):
    #Hardcoded at N=418
    N = 418
    coords = np.empty([418,2])
    dex = 0
    for i in range(22):
        rise = np.sqrt(3)/2*i
        for j in range(19):
            if i%2 == 0: offset = 0.5
            else: offset = 0.0
            coords[dex] = [offset+j,rise]
            dex += 1
    #print(coords) 
    
    #To prevent some mild annoyances with the wrap function
    coords += 0.001

    if a != 0:
        for i in range(coords.shape[0]):
            coords[i] += (np.random.rand(2) - 0.5) * a

    latmat = np.array([[19,0],[0,19]])
    wrap(coords, latmat)

    return(coords, latmat)


## Load (Ge)
#This assumes that we know the particle number
def gen_Ge(path):
    with open(path, 'r') as f:
        f.readline() #Tossing line that contains the dimension
        lv1 = f.readline().split(' ')[:2]
        lv2 = f.readline().split(' ')[:2]
        coords = []
        for line in f:
            
            coords.append(line.split(' ')[:2])

    coords = np.array(coords).astype(float)
    latmat = np.array([lv1,lv2]).T.astype(float)
    return(coords, latmat)

## Load (TJ)
def gen_TJ(path):
    with open(path, 'r') as f:
        TJtype = f.readline().split('\t')[-1]
        f.readline()
        f.readline()
        if TJtype == 'mono\n': f.readline()
        
        lv1 = f.readline().split('\t')[:2]
        lv2 = f.readline().split('\t')[:2]
        f.readline()

        coords = []
        rads = []

        for line in f:
            tline = line.split('\t')
            coords.append(tline[:2])
            rads.append(0)
            #rads.append(tline[2])            

    coords = np.array(coords).astype(float)
    rads = np.array(rads).astype(float)
    latmat = np.array([lv1,lv2]).astype(float).T

    return(coords, latmat, rads)

###Misc
## Make supercell
def tile(coords, latmat):
    tcoords = np.zeros([coords.shape[0]*9, 2])
    imct = 0
    for v1 in [0,-1,1]:
        for v2 in [0,-1,1]:
            im = np.copy(coords)+v1*latmat[:,0].T+v2*latmat[:,1].T 
            tcoords[imct*coords.shape[0]:(imct+1)*coords.shape[0]] = im
            imct += 1
    return(tcoords)

def wrap(coords, latmat):
    for i in range(coords.shape[0]):
        lcoord = np.matmul(np.linalg.inv(latmat),coords[i])
        if lcoord[0] < 0: lcoord[0] += 1
        if lcoord[0] > 1: lcoord[0] -= 1
        if lcoord[1] < 0: lcoord[1] += 1
        if lcoord[1] > 1: lcoord[1] -= 1
        coords[i] = np.matmul(latmat,lcoord)
        
###################################################################################################
### Functions to make tessellations (2D)

## Various Tilings
def tile_gen(name):
    if 'Z2D' in name:
        #Hard coding 20x20

        defects = float(name.split('_')[-1])
        defects = int(defects * 19*19)

        dlocs = np.random.choice(19*19, size=defects, replace=False)

        coords, toss = gen_URL(400,0.0)

        edges = []
        for i in range(20):
            for j in range(20):
                if i != 19: edges.append([i*20+j,(i+1)*20+j])
                if j != 19: edges.append([i*20+j,i*20+j+1])

        for d in dlocs:
            boxloc = (d//19)*20 + d%19
            if random.choice([True, False]):
                edges.append([boxloc,boxloc+21])
            else: edges.append([boxloc+20,boxloc+1])

        return(coords, np.array(edges))

    if name == 'Hex':
        coords = np.zeros([420,2])
        for j in range(21):
            coords[j] = [(0.5*((j+1)%2)),np.sqrt(3)/2*j]
            coords[21+j] = [(0.5*((j)%2))+1.5,np.sqrt(3)/2*j]
        

        for i in range(1,10):
            coords[42*i:42*(i+1)] = coords[42*(i-1):42*(i)] + [3,0]

        edges = []
        for i in range(20):
            for j in range(21):
                if i%2 == 0 and i != 20-1:
                    if j%2 == 0: edges.append([i*21+j,(i+1)*21+j])
                elif i != 20-1:
                    if (j+1)%2 == 0: edges.append([i*21+j,(i+1)*21+j])
                if j != 20: edges.append([i*21+j,i*21+j+1])
        return(coords, np.array(edges))

    if name == 'Kago':
        coords = np.zeros([384,2])
        for i in range(16):
            coords[i] = [0.5+i,0]
        for i in range(9):
            coords[i+16] = [i*2,np.sqrt(3)/2]
        coords[16+9:16+16+9] = coords[:16]+[0,np.sqrt(3)]
        for i in range(8):
            coords[i+41] = [1+(2*i),3/2*np.sqrt(3)]
        
        
        for i in range(1,7):
            coords[i*49:(i+1)*49] = coords[:49]+[0,2*np.sqrt(3)*i]
        
        coords[-41:] = coords[:41]+[0,14*np.sqrt(3)]

        plt.scatter(coords[:,0],coords[:,1])
        plt.gca().set_aspect('equal')
        #plt.show()

        edges = []
        
        for i in range(15):
            edges.append([i,i+1])
        edges.append([16,0])
        edges.append([16,25])
        for i in range(1,8):
            edges.append([16+i,2*i-1])
            edges.append([16+i,(2*i-1)+25])
            edges.append([16+i,2*i])
            edges.append([16+i,(2*i)+25])
        edges.append([24,15])
        edges.append([24,40])
        for i in range(25,40):
            edges.append([i,i+1])
        for i in range(8):
            edges.append([41+i,(2*i+1)+25])
            edges.append([41+i,(2*i+1)+49])
            edges.append([41+i,(2*i)+25])
            edges.append([41+i,(2*i)+49])
        
        edge_up = []
        for t in range(7):
            for i in edges:
                edge_up.append([i[0]+t*49, i[1]+t*49])
        for i in edges[:62]:
            edge_up.append([i[0]+7*49, i[1]+7*49])

        return(coords, np.array(edge_up))


## Delaunay
def delaunay_graph_2d(points):
    # Perform Delaunay triangulation
    delaunay = sp.spatial.Delaunay(points)
    triangles = np.array(delaunay.simplices)

    # Convert triangles to edges
    edges = restructure_array_numpy(triangles)

    # Remove duplicate edges
    unique_edges = remove_duplicate_rows(edges)

    return unique_edges

def restructure_array_numpy(triangles):
    reshaped = triangles
    output = np.empty((len(reshaped) * 3, 2), dtype=triangles.dtype)
    output[::3, 0] = reshaped[:, 0]
    output[::3, 1] = reshaped[:, 1]
    output[1::3, 0] = reshaped[:, 1]
    output[1::3, 1] = reshaped[:, 2]
    output[2::3, 0] = reshaped[:, 2]
    output[2::3, 1] = reshaped[:, 0]
    return output

def remove_duplicate_rows(arr):
    dtype = [('min', arr.dtype), ('max', arr.dtype)]
    structured = np.empty(arr.shape[0], dtype=dtype)
    structured['min'] = np.minimum(arr[:, 0], arr[:, 1])
    structured['max'] = np.maximum(arr[:, 0], arr[:, 1])
    _, unique_indices = np.unique(structured, return_index=True)
    unique_indices.sort()
    return arr[unique_indices]

def gabriel_graph_2d(points):
    """
    Compute the Gabriel graph of a set of points in 2D or 3D space.

    Parameters:
    - points: (N, D) numpy array where N is the number of points and D is the dimension (2 or 3).

    Returns:
    - edge_array: (M, 3) numpy array where each row represents an edge in the format [node1, node2, weight].
    """

    # Build a KD-tree for efficient neighbor searches
    tree = sp.spatial.cKDTree(points)

    # Compute the Delaunay triangulation
    tri = sp.spatial.Delaunay(points)
    simplices = tri.simplices  # Indices of points forming the simplices


    # Generate all possible edges from the simplices
    if simplices.shape[1] == 3:  # 2D case
        edges = np.vstack([simplices[:, [0, 1]],
                           simplices[:, [1, 2]],
                           simplices[:, [2, 0]]])

    elif simplices.shape[1] == 4:  # 3D case
        edges = np.vstack([simplices[:, [0, 1]],
                           simplices[:, [0, 2]],
                           simplices[:, [0, 3]],
                           simplices[:, [1, 2]],
                           simplices[:, [1, 3]],
                           simplices[:, [2, 3]]])
    else:
        raise ValueError('Input points must be 2D or 3D.')



    # Sort and remove duplicate edges
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)


    # Compute midpoints and radii for the Gabriel condition
    i = edges[:, 0]
    j = edges[:, 1]
    points_i = points[i]
    points_j = points[j]
    midpoints = (points_i + points_j) / 2
    radii = 0.5 * np.linalg.norm(points_i - points_j, axis=1)
    weights = 2 * radii  # Edge weights are the Euclidean distances

    # Query KD-tree to find neighboring points within the radius
    idx_list = tree.query_ball_point(midpoints, radii)

    # Build the edge list for the Gabriel graph
    edge_list = []
    for k in range(len(edges)):
        idx = set(idx_list[k]) - {i[k], j[k]}
        if not idx:
            edge_list.append([i[k], j[k], weights[k]])

    edge_array = np.array(edge_list)[:,:2]
    return edge_array.astype(int)


def process_edges_and_points(points, edges):
    # Step 1: Sort points by first column
    sorted_points, sorted_indices = sort_points_by_first_column(points)

    # Step 2: Compute centroids and radii (half-lengths) of edges
    centroids, radii = compute_edge_centroids_and_half_lengths(edges, points)

    # Step 3: Check which edges have no points within their radii
    new_edge_indices = check_points_in_centroids(sorted_points, centroids, radii)

    return new_edge_indices

def compute_edge_centroids_and_half_lengths(edges, points):
    start_points = points[edges[:, 0]]
    end_points = points[edges[:, 1]]

    centroids = (start_points + end_points) / 2
    diff_vectors = end_points - start_points
    half_lengths = np.sqrt(np.sum(diff_vectors**2, axis=1)) / 2

    return centroids, half_lengths


def sort_points_by_first_column(points):
    sorted_indices = np.argsort(points[:, 0])
    sorted_points = points[sorted_indices]
    return sorted_points, sorted_indices

def check_points_in_centroids(sorted_points, centroids, radii):
    new_edges = []
    n_centroids = len(centroids)

    for i in range(n_centroids):
        centroid = centroids[i]
        radius = radii[i]

        left_index = np.searchsorted(sorted_points[:, 0], centroid[0] - radius)
        right_index = np.searchsorted(sorted_points[:, 0], centroid[0] + radius, side='right')

        potential_points = sorted_points[left_index:right_index]

        y_mask = np.abs(potential_points[:, 1] - centroid[1]) < radius
        close_points = potential_points[y_mask]

        if len(close_points) == 0:
            new_edges.append(i)
            continue

        distances = np.sqrt(np.sum((close_points - centroid)**2, axis=1))

        if not np.any(distances < radius):
            new_edges.append(int(i))

    return np.array(new_edges)

def delaunay_to_unique_edges(points):
    # Perform Delaunay triangulation
    delaunay = sp.spatial.Delaunay(points)
    triangles = np.array(delaunay.simplices)

    # Convert triangles to edges
    edges = restructure_array_numpy(triangles)

    # Remove duplicate edges
    unique_edges = remove_duplicate_rows(edges)

    return unique_edges

## Voronoi
#Omitting this for the time being
def voronoi_graph_2d(points):
    vor = sp.spatial.Voronoi(points)
    using = vor.ridge_vertices
    using = remove_duplicate_rows(np.array(using))
    using = using[np.argwhere(np.sum(using>=0,axis=1)==2).flatten()]
    newpoints = vor.vertices
    return(newpoints, np.array(using))

## Delaunay-Centroidal
def delcent_graph_2d(points):
    tri = sp.spatial.Delaunay(points)
    newpoints = np.zeros([tri.simplices.shape[0],2])
    for i in range(newpoints.shape[0]):
        newpoints[i] = np.mean(points[tri.simplices[i]],axis = 0)
    edges = []
    for i in range(tri.simplices.shape[0]):
        for j in range(3):
            if tri.neighbors[i][j] != -1: edges.append([i,tri.neighbors[i][j]])
    return(newpoints, np.array(edges))

###################################################################################################
### Functions to make the graphs with various boundary conditions
def clip_boundary(alledge,tcoords,latmat):
    outvec = bdry_check(tcoords,latmat) 
    #First, check if an edge straddles the box boundary
    edp = modedge_find(alledge, outvec, tcoords)
    
    segs = np.array([[[0,0],[latmat[0,0],latmat[1,0]]],
            [[latmat[0,0],latmat[1,0]],[latmat[0,1]+latmat[0,0],latmat[1,1]+latmat[1,0]]],
            [[latmat[0,1],latmat[1,1]],[latmat[0,1]+latmat[0,0],latmat[1,1]+latmat[1,0]]],
            [[0,0],[latmat[0,1],latmat[1,1]]]])
    

    dex = len(edp) - 1

    enddex = tcoords.shape[0]
    newedges = []
    for i in edp[::-1]:
        if len(i) == 4:
            newedge = [i[0],i[1],0]
            fixdex = i[-1]
            clippoint = [None, None] 
            segdex = 0
            while clippoint[0] == None:
                clippoint = intersect(tcoords[i[0]],tcoords[i[1]],segs[segdex,0],segs[segdex,1])
                segdex += 1
            clippoint = np.array([clippoint])
            tcoords = np.append(tcoords,clippoint,axis=0)
            
            ###No, need to add in new points, look at how I did this for Chenxi...
            newedge[fixdex] = tcoords.shape[0] - 1
            #tcoords[i[fixdex]] = clippoint
            #havefixed.append(i[fixdex])
            newedge[2] = np.linalg.norm(tcoords[i[0]]-tcoords[i[1]])
            edp.pop(dex)
            newedges.append(newedge)
            #edp[dex] = edp[dex][:3]
        dex -= 1
    edp = np.array(edp)
    edp = remove_duplicate_rows(edp)
    edp = np.concatenate((edp,np.array(newedges)),axis=0)
    vertindexes = np.unique(edp[:,:2])
    fixedge = np.copy(edp)
    for i in range(fixedge.shape[0]):
        fixedge[i,0] = np.argwhere(vertindexes==edp[i,0])
        fixedge[i,1] = np.argwhere(vertindexes==edp[i,1])

    vertlocs = np.zeros([vertindexes.shape[0],2])
    for i in range(vertlocs.shape[0]):
        vertlocs[i]=tcoords[int(vertindexes[i])]
    return(vertlocs, fixedge)

def delete_boundary(alledge,tcoords,latmat):
    outvec = bdry_check(tcoords,latmat)
    #First, check if an edge straddles the box boundary
    edp = modedge_find(alledge, outvec, tcoords) 
    
    dex = len(edp)-1
    for i in edp[::-1]:
        if len(i) == 4:
            edp.pop(dex)     
        dex -= 1

    edp = np.array(edp)
    edp = remove_duplicate_rows(edp)
    vertindexes = np.unique(edp[:,:2])
    fixedge = np.copy(edp)
    for i in range(fixedge.shape[0]):
        fixedge[i,0] = np.argwhere(vertindexes==edp[i,0])
        fixedge[i,1] = np.argwhere(vertindexes==edp[i,1])

    vertlocs = np.zeros([vertindexes.shape[0],2])
    for i in range(vertlocs.shape[0]):
        vertlocs[i]=tcoords[int(vertindexes[i])]
    return(vertlocs, fixedge)


def prune_boundary(alledge,tcoords,latmat):
    outvec = bdry_check(tcoords,latmat)
    #First, check if an edge straddles the box boundary
    edp = modedge_find(alledge, outvec, tcoords)

    dex = len(edp)-1
    for i in edp[::-1]:
        if len(i) == 4:
            edp.pop(dex)
        dex -= 1

    edp = np.array(edp)
    edp = remove_duplicate_rows(edp)

    edgedexs, dexcounts = np.unique(edp[:,:2].flatten(),return_counts=True)
    deldex = np.argwhere(dexcounts == 1)
    deg1 = edgedexs[deldex].astype(int)

    while deg1.shape[0] > 0:
        for todel in deg1:
            tdex = np.argwhere(edp == todel)[0,0]
            edp = np.delete(edp,tdex,axis=0)
        edgedexs, dexcounts = np.unique(edp[:,:2].flatten(),return_counts=True)
        deldex = np.argwhere(dexcounts == 1)
        deg1 = edgedexs[deldex].astype(int)

    vertindexes = np.unique(edp[:,:2])
    fixedge = np.copy(edp)
    for i in range(fixedge.shape[0]):
        fixedge[i,0] = np.argwhere(vertindexes==edp[i,0])
        fixedge[i,1] = np.argwhere(vertindexes==edp[i,1])


    vertlocs = np.zeros([vertindexes.shape[0],2])
    for i in range(vertlocs.shape[0]):
        vertlocs[i]=tcoords[int(vertindexes[i])]
    return(vertlocs, fixedge)

def bdry_check(tcoords,latmat):
    #Find vertices that lie outside the fundamental cell
    lcoords = np.copy(tcoords)
    latinv = np.linalg.inv(latmat)
    for i in range(lcoords.shape[0]):
        lcoords[i] = np.matmul(latinv,lcoords[i])
    lcoords = np.abs(np.floor(lcoords))
    lsum = np.sum(lcoords,axis=1).astype('bool') 
    
    return(lsum)

def modedge_find(alledge, outvec, tcoords):
    trimmededge = []

    for e in alledge:
        length = np.linalg.norm(tcoords[e[0]]-tcoords[e[1]])
        if outvec[e[0]] == False and outvec[e[1]] == False:
            trimmededge.append([e[0],e[1],length])
        elif outvec[e[0]] == True and outvec[e[1]] == True:
            continue
        else:
            if outvec[e[0]] == True: trimmededge.append([e[0],e[1],length,0]) 
            else: trimmededge.append([e[0],e[1],length,1])

    return(trimmededge)

def intersect(p1, p2, p3, p4):
    """Checks if line segment p1-p2 intersects with line segment p3-p4"""
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    if denom == 0:  # parallel lines
        return [None,None]

    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denom
    if ua < 0 or ua > 1:  # intersection point not on line segment p1-p2
        return [None,None]

    ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom
    if ub < 0 or ub > 1:  # intersection point not on line segment p3-p4
        return [None,None]

    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)
    return np.array([x, y])

###################################################################################################
### Functions for visualization
def configvis(coords, latmat):
    plt.scatter(coords[:,0], coords[:,1])
    plt.xlim(np.min(latmat[0,:]),np.max(latmat[0,:]))
    plt.ylim(np.min(latmat[1,:]),np.max(latmat[1,:]))
    plt.show()

def netvis(edges, coords, latmat,color):
    for e in edges:
        e = e[:2].astype(int)
        plt.plot(coords[e,0],coords[e,1], color)
    plt.plot([0,latmat[0,0],latmat[0,0]+latmat[0,1],latmat[0,1],0],[0,latmat[1,0],latmat[1,0]+latmat[1,1],latmat[1,1],0],'k')
    plt.xlim(np.min(latmat[0,:]*-1),np.max(latmat[0,:]*2))
    plt.ylim(np.min(latmat[1,:]*-1),np.max(latmat[1,:]*2))
    plt.gca().set_aspect('equal')
    #plt.show()

###################################################################################################
### Script 
conftype = str(sys.argv[1])
N = int(sys.argv[2])
tess = str(sys.argv[3])
bdry = str(sys.argv[4])
name = str(sys.argv[5])

if conftype == 'URL':
    ct = 'URL'
    a = float(sys.argv[6])
    params = 'a' + str(a)
if conftype == 'A2':
    ct = 'A2'
    a = float(sys.argv[6])
    params = 'a' + str(a)
if conftype == 'Ge' or conftype == 'TJ':
    ct = conftype 
    path = sys.argv[6]
    params = sys.argv[7]

if conftype == 'URL':
    coords, latmat = gen_URL(N,a) 
elif conftype == 'A2':
    coords, latmat= gen_A2(N,a)
elif conftype == 'Ge':
    coords, latmat = gen_Ge(path)
elif conftype == 'TJ':
    coords, latmat, rads = gen_TJ(path)
elif conftype == 'Ti':
    ct = 'Ti'
    ttype = str(sys.argv[6])
    params = ttype
    coords, alledge = tile_gen(ttype)


if conftype != 'Ti': tcoords = tile(coords,latmat)

if tess == 'D':
    alledge = delaunay_graph_2d(tcoords)
elif tess == 'V':
    tcoords, alledge = voronoi_graph_2d(tcoords)
elif tess == 'C':
    tcoords, alledge = delcent_graph_2d(tcoords)
elif tess == 'G':
    alledge = gabriel_graph_2d(tcoords)

if bdry == 'delete':
    vertlocs, fixedge = delete_boundary(alledge, tcoords, latmat)
elif bdry == 'prune':
    vertlocs, fixedge = prune_boundary(alledge, tcoords, latmat)
elif bdry == 'clip':
    vertlocs, fixedge = clip_boundary(alledge, tcoords, latmat)
else:
    vertlocs = coords
    fixedge = np.hstack((alledge,np.ones((alledge.shape[0],1))))
    for i in range(fixedge.shape[0]):
        fixedge[i,-1] = np.linalg.norm(vertlocs[int(fixedge[i,0])]-vertlocs[int(fixedge[i,1])])


filetag = ct + "_" + params + "_" + tess + "_N" + str(N) + "_" + bdry + "_" + name

if bdry != 'periodic':
    padverts = np.hstack((vertlocs,np.zeros((vertlocs.shape[0],1))))
    np.savetxt("./configs/configs/vertlocs_"+filetag+".txt",padverts)
    np.savetxt("./configs/edgelists/edges_"+filetag+".txt",fixedge)
else:
    torusverts = FC_to_torus(vertlocs,latmat)
    np.savetxt("./configs/configs/vertlocs_"+filetag+".txt",torusverts)
    np.savetxt("./configs/edgelists/edges_"+filetag+".txt",fixedge)




