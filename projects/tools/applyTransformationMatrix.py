import numpy as np
import msh
import sys

def readMatrixFile(path):
    MAT = np.zeros(shape=(4,4))
    with open(path) as f:
        for i,l in enumerate(f.readlines()):
            if i<4:
                elts = [float(x) for x in l.strip().split()]
                MAT[i] = elts
    print MAT
    return np.array(MAT)

if __name__=="__main__":
    mesh = msh.Mesh(sys.argv[1])
    MAT = readMatrixFile(sys.argv[2])
    mesh.verts = np.array([ np.insert(  np.dot( MAT, np.append(v[:3],[1]) )[:3], 3, v[-1]  ) for v in mesh.verts ])
    mesh.write("out.mesh")
