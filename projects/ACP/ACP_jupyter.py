import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import random
import math

import icp
import pca
import msh

def writeDiff(d1, d2, tris, fileroot):
    with open(fileroot+".mesh","w") as f:
        f.write("MeshVersionFormatted 2\nDimension 3\nVertices\n"+str(2*len(d1))+"\n")
        for i,d in enumerate(d1):
            f.write(" ".join([str(x) for x in d]) + " 1\n")
        f.write("Triangles\n"+str(len(tris))+"\n")
        for t in tris:
            f.write(" ".join([str(x+1) for x in t]) + " 0\n")
    with open(fileroot+".sol","w") as f:
        f.write("MeshVersionFormatted 1\nDimension 3\nSolAtVertices\n"+str(len(d1))+"\n1 1\n")
        for i,d in enumerate(d1):
            f.write(str(np.linalg.norm(d - d2[i])) + "\n")
def writeDefault(d1,d2):
    with open("DEFAULT.medit","w") as f:
        f.write("WindowSize 400 600\n")
        f.write("Palette ")
        dmax = np.max([np.linalg.norm(x-y) for x,y in zip(d1,d2)])
        f.write(" ".join( [str(i*dmax/4.0) for i in range(5)] ))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
def rotate_verts(verts, angleRandomness=0.3, tranlsationRandomness=1):
    axis = [random.random(),random.random(),random.random()]
    theta = angleRandomness*random.random()
    tr = [tranlsationRandomness*random.random() for i in range(3)]
    newV = np.array([np.dot(rotation_matrix(axis,theta), v) for v in verts])+tr
    return newV
#v = [3, 5, 0]
#axis = [4, 4, 1]
#theta = 1.2
#print(np.dot(rotation_matrix(axis,theta), v))

if __name__ == "__main__":

    verts,tris = msh.readMesh("/Users/guestaccount/dev/ACP/templateSurf.mesh")
    #Avec les donnees brutes
    directory = "data/"
    exclude   = ["MassDisp.21.sol.txt","MassDisp.24.sol.txt","MassDisp.48.sol.txt","MassDisp.54.sol.txt"]
    files     = [[ int(f.split(".")[1]), directory+f ] for f in os.listdir(directory) if ".sol.txt" in f and f not in exclude]
    data      = [pca.read_data(f[1]) for f in sorted(files)]
    #alignedData = icp.align_data(data,data[-1],offset=verts, write=True)
    nb        = int(sys.argv[1]) if len(sys.argv)>1 else 1
    X         = pca.PCA(data[:-1], data[-1], nb)
    writeDiff(X,data[-1],tris,"rec")

    XX = X-data[-1]
    #os.system("cp templateSurf.mesh results/reconstruction.mesh")
    #msh.write_solution(X, "results/reconstruction.sol")
    #msh.writeMesh(verts+X,tris,"results/reconstruction.mesh")
    #msh.write_solution([[np.linalg.norm(d-r)] for d,r in zip(data[-1],X)], "results/reconstruction.sol")

    directory = "data_aligned/"
    exclude   = []
    files     = [[ int(f.split(".")[1]), directory+f ] for f in os.listdir(directory) if ".sol.txt" in f and f not in exclude]
    data      = [pca.center_data(pca.read_data(f[1])) for f in sorted(files)]
    nb        = int(sys.argv[1]) if len(sys.argv)>1 else 1
    Y         = pca.PCA(data[:-1], data[-1], nb)
    writeDiff(verts+Y,verts+data[-1],tris,"ali")
    YY = Y-data[-1]

    writeDefault(XX,YY)
    #os.system("cp templateSurf.mesh results/alignement.mesh")
    #msh.write_solution(Y, "results/alignement.sol")
    #msh.writeMesh(verts+Y,tris,"results/alignement.mesh")
