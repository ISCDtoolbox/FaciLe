from scipy   import ndimage as nd
from skimage import measure as mea
import numpy as np
import os
import msh
import sys

"""
os.system("powercrust -i " + root + ".xyz > log.txt 2>&1")
os.system("LC_ALL=C meshlab pc.off >> log.txt 2>&1")
os.system("boundingmesh -v 200 -d Outward pc.off >> log.txt 2>&1")
"""

def readObj(inFile):
    verts, tris = [],[]
    with open(inFile) as f:
        for l in f.readlines():
            L = l.strip()
            if len(L):
                L = L.split()
                if L[0] == "v":
                    verts.append([float(x) for x in L[1:]] + [0])
                if L[0] == "f":
                    tris.append([int(x)-1 for x in L[1:]] + [0])
    return np.array(verts), np.array(tris)
def readOff(offFile):
    with open(offFile) as f:
        LINES = f.readlines()
        if 'OFF' not in LINES[0]:
            print 'Not a valid OFF header'
        n_verts, n_faces, n_dontknow = tuple([int(s) for s in LINES[1].strip().split(' ')])
        verts = []
        for i_vert in range(n_verts):
            verts.append([float(s) for s in LINES[i_vert+2].strip().split()])
        faces = []
        for i_face in range(n_faces):
            faces.append([int(s) for s in LINES[n_verts + i_face+2].strip().split()[1:]])
        return np.array(verts), np.array(faces)
def binaryToXYZ(data, xyzFile):
    X,Y,Z = data.nonzero()
    with open(xyzFile,"w") as f:
        for x,y,z in zip(X,Y,Z):
            f.write(" ".join([str(x), str(y), str(z)]) + "\n")

def coorToInd(pt,res,xmin,xmax,ymin,ymax,zmin,zmax):
    pt[0] = 5+int( res*(pt[0]-xmin)/(xmax-xmin) - 1e-8)
    pt[1] = 5+int( res*(pt[1]-ymin)/(ymax-ymin) - 1e-8)
    pt[2] = 5+int( res*(pt[2]-zmin)/(zmax-zmin) - 1e-8)
    return pt
def addCube(xmin,xmax,ymin,ymax,zmin,zmax,ref=100):
    cubeVerts = np.array([
        [xmin, ymin, zmin],
        [xmin, ymin, zmax],
        [xmax, ymin, zmin],
        [xmax, ymin, zmax],
        [xmin, ymax, zmin],
        [xmin, ymax, zmax],
        [xmax, ymax, zmin],
        [xmax, ymax, zmax]
    ])
    cubeTris = np.array([
        [0,1,2],
        [1,3,2],
        [4,6,5],
        [5,6,7],
        [1,5,3],
        [3,5,7],
        [2,6,4],
        [0,2,4],
        [3,7,6],
        [2,3,6],
        [0,4,1],
        [1,4,5]
    ])
    return np.insert(cubeVerts,3,ref,axis=1), np.insert(cubeTris,3,ref,axis=1)
def ptsToXYZCubes(pts, res=100):
    data = np.zeros(shape=(res+10,res+10,res+10), dtype="bool_")
    xmin,xmax = np.min(pts[:,0]), np.max(pts[:,0])
    ymin,ymax = np.min(pts[:,1]), np.max(pts[:,1])
    zmin,zmax = np.min(pts[:,2]), np.max(pts[:,2])
    resultingScale = [(xmax-xmin)/res, (ymax-ymin)/res, (zmax-zmin)/res]
    tmpPts = np.array([coorToInd(pt,res,xmin,xmax,ymin,ymax,zmin,zmax) for pt in np.copy(pts)])
    tmpPts = tmpPts.astype(np.int16, copy=False)
    for pt in tmpPts:
        data[pt[0], pt[1], pt[2]] = True
    return data, resultingScale

def indInfCorner(i,res):
    return [res-1+j-i for j in range(i+1)], [j for j in range(i+1)]
def indInfCenter(i,res):
    return [i+j+1 for j in range(res-1-i)], [j for j in range(res-1-i)]
def indSupCorner(i,res):
    return [j for j in range(i+1)], [res-1-i+j for j in range(i+1)]
def indSupCenter(i,res):
    return [j for j in range(res-1-i)], [i+j+1 for j in range(res-1-i)]
def carveDiag(data, newData, tr=(0,1,2), inv=False):
    #En diagonale
    #Diagonales en 2d:
    #Tableau de taille res + 1, avec chaque item de longueur 5
    #m[0]:diagonale du coin inferieur+de la bande a cote de la diagonale
    res = len(data)
    R = (res-1)/2

    if not inv:

        iInfCor = [indInfCorner(i,res) for i in range(R)]
        iInfCen = [indInfCenter(i,res) for i in range(R)]
        iSupCor = [indSupCorner(i,res) for i in range(R)]
        iSupCen = [indSupCenter(i,res) for i in range(R)]
        iDiag   = [range(res), range(res)]


        for d,n in zip(np.transpose(data,tr), np.transpose(newData,tr)):
            for i in range(R):
                has = np.any(d[iInfCor[i]])
                if not has:
                    n[iInfCor[i]] = has

                has = np.any(d[iInfCen[i]])
                if not has:
                    n[iInfCen[i]] = has
            for i in range(R):
                has = np.any(d[iSupCor[i]])
                if not has:
                    n[iSupCor[i]] = has

                has = np.any(d[iSupCen[i]])
                if not has:
                    n[iSupCen[i]] = has
            has = np.any(d[iDiag])
            if not has:
                n[iDiag] = has
    else:
        iInfCor = [indInfCorner(i,res) for i in range(R)]
        iInfCen = [indInfCenter(i,res) for i in range(R)]
        iSupCor = [indSupCorner(i,res) for i in range(R)]
        iSupCen = [indSupCenter(i,res) for i in range(R)]
        iDiag   = [range(res), range(res)]

        iInfCor = [[i[0],res-1-np.array(i[1])] for i in iInfCor ]
        iInfCen = [[i[0],res-1-np.array(i[1])] for i in iInfCen ]
        iSupCor = [[i[0],res-1-np.array(i[1])] for i in iSupCor ]
        iSupCen = [[i[0],res-1-np.array(i[1])] for i in iSupCen ]
        iDiag   = [range(res), range(res)[::-1] ]

        for d,n in zip(np.transpose(data,tr), np.transpose(newData,tr)):
            for i in range(R):
                has = np.any(d[iInfCor[i]])
                if not has:
                    n[iInfCor[i]] = has

                has = np.any(d[iInfCen[i]])
                if not has:
                    n[iInfCen[i]] = has
            for i in range(R):
                has = np.any(d[iSupCor[i]])
                if not has:
                    n[iSupCor[i]] = has

                has = np.any(d[iSupCen[i]])
                if not has:
                    n[iSupCen[i]] = has
            has = np.any(d[iDiag])
            if not has:
                n[iDiag] = has

    return newData
def carveAxis(data, newData,tr=(0,1,2)):
    mask = np.sum(np.transpose(data,tr),axis=0)>0
    mask = nd.binary_closing(mask,iterations=3)
    for d in np.transpose(newData, tr):
        d[mask==False] = False
    return newData
def get3DdiagonalIndices(res):
    #First half and middle slice
    inds = []
    for r in range(res):
        for i in range(r+1):
            for j in np.arange(0,r+1-i):
                for k in np.arange(r-i-j,r+1-i-j):
                    if i==0 or j==0 or k==0:
                        inds.append([i,j,k])
    #Get the triangles
    for i in range(res-1):
        for j in range(res-1):
            if i+j<= res-2:
                inds.append([0,res-1-i,res-1-j])
                inds.append([res-1-i,0,res-1-j])
                inds.append([res-1-i,res-1-j,0])
    inds = np.array(inds)
    #Creating the diagonals
    offset = np.max(inds,axis=1)
    diags  = [ np.array([i+j for j in range(res - o)]) for i,o in zip(inds,offset) ]
    return diags
def carveSuperDiagonal(data, newData):
    diags = get3DdiagonalIndices(len(data))
    for d in diags:

        transformedData = data
        transformedNewData = newData
        if not np.any(transformedData[d[:,0],d[:,1],d[:,2]]):
            transformedNewData[d[:,0],d[:,1],d[:,2]]=False

        transformedData = data[::-1,:,:]
        transformedNewData = newData[::-1,:,:]
        if not np.any(transformedData[d[:,0],d[:,1],d[:,2]]):
            transformedNewData[d[:,0],d[:,1],d[:,2]]=False

        transformedData = data[:,::-1,:]
        transformedNewData = newData[:,::-1,:]
        if not np.any(transformedData[d[:,0],d[:,1],d[:,2]]):
            transformedNewData[d[:,0],d[:,1],d[:,2]]=False

        transformedData = data[:,:,::-1]
        transformedNewData = newData[:,:,::-1]
        if not np.any(transformedData[d[:,0],d[:,1],d[:,2]]):
            transformedNewData[d[:,0],d[:,1],d[:,2]]=False

def spaceCarve(data):
    newData = np.invert(np.zeros(shape=data.shape, dtype="bool_"))

    #Super Diagonal
    carveSuperDiagonal(data, newData)

    #Diagonal axis
    carveDiag(data, newData, (0,1,2))
    carveDiag(data, newData, (0,2,1), inv=True)
    carveDiag(data, newData, (1,0,2))
    carveDiag(data, newData, (1,2,0), inv=True)
    carveDiag(data, newData, (2,0,1))
    carveDiag(data, newData, (2,1,0), inv=True)

    #Principal axis
    carveAxis(data, newData, (0,1,2))
    carveAxis(data, newData, (1,0,2))
    carveAxis(data, newData, (2,1,0))

    return newData

def command(cmd, displayOutput=False):
    err = os.system(cmd) if displayOutput else os.system(cmd + " > tmp_out.txt 2>tmp_err.txt")
    if err:
        print "An error happened while executing:\n"+cmd+"\nLook in tmp_out.txt or tmp_err.txt"
        sys.exit()
    else:
        os.system("rm tmp_out.txt tmp_err.txt >/dev/null 2>&1")

resolution = 71

if __name__=="__main__":
    #Parsing arguments
    root = sys.argv[1].split("/")[-1][:-5]


    print "1 - Converting .mesh to binary .xyz"
    print "-  1.1 - Opening the mesh file"
    mesh = msh.Mesh(sys.argv[1])
    print "-  1.2 - Converting to binary point data"
    binaryData, totalScale = ptsToXYZCubes(mesh.verts,resolution)


    print "2 - Creating the filled volume"
    print "-  2.1 - Space carving"
    newData = spaceCarve(binaryData)
    newData = nd.binary_closing(newData,structure=nd.generate_binary_structure(3, 3),iterations=3)
    print "-  2.2 - Marching cubes"
    verts, faces = mea.marching_cubes(volume=newData, level=0.5)
    recon = msh.Mesh()
    recon.verts = np.insert(np.array(verts),3,0,axis=1)
    recon.tris = np.insert(np.array(faces),3,0,axis=1)
    recon.computeBBox()
    print "-  2.3 - Writing the scaled reconstructed .mesh file"
    recon.fitTo(mesh,keepRatio=False)
    recon.write(root + ".hull.mesh")
    del binaryData, newData, recon
    print "2.4 - Remeshing the hull"
    os.system("mmgs_O3 " + root + ".hull.mesh -nr -hausd " + str(np.max(mesh.dims)/100))
