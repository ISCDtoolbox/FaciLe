from scipy   import ndimage as nd
from skimage import measure as mea
import numpy as np
import os
import msh
import sys
import argparse

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

def parse():
    parser = argparse.ArgumentParser(description="Creates an adequate shell for warping")
    parser.add_argument("-i", "--input", help="input .mesh file", type=str, required=True)
    parser.add_argument("-o", "--output", help="transformed .mesh file", type=str, required=True)
    parser.add_argument("-r", "--resolution", help="resolution of the carving algorithm, has to be odd", type=int, default=51)
    parser.add_argument("-m", "--remesh", help="Apply a final surface remeshing to the carved mesh", action="store_true")
    return parser.parse_args()

def checkArgs(args):
    if not os.path.isfile(args.input):
        print args.input + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.input)[1] == ".mesh":
        print args.input + " is not a .mesh file"
        sys.exit()
    if not os.path.splitext(args.output)[1] == ".mesh":
        print "Output file must be in the .mesh format"
        sys.exit()
    if args.resolution < 11 or args.resolution>301:
        print "The resolution must be in [11, 301]"
        sys.exit()

if __name__=="__main__":
    args = parse()
    checkArgs(args)

    print "1 - Converting .mesh to binary .xyz"
    print "-  1.1 - Opening the mesh file"
    mesh = msh.Mesh(args.input)
    print "-  1.2 - Converting to binary point data"
    binaryData, totalScale = ptsToXYZCubes(mesh.verts,args.resolution)


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
    recon.write(args.output)
    del binaryData, newData, recon
    if args.remesh:
        print "2.4 - Remeshing the hull"
        os.system("mmgs_O3 " + args.output + " -nr -hausd " + str(np.max(mesh.dims)/100) + " -o " + args.output)
