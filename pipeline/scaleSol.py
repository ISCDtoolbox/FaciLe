import sys
sys.path.append("../projects/tools")
import msh
import os
import numpy as np
from scipy.spatial.distance import cdist

def nearest_neighbor(src, dst):
    all_dists = cdist(src, dst, 'euclidean')
    indices = all_dists.argmin(axis=1)
    distances = all_dists[np.arange(all_dists.shape[0]), indices]
    return distances, indices

if __name__=="__main__":
    one = False
    two = True

    if one:
        directory = "/Data/Facile2/muscles/"
        files = [ os.path.join(directory, f) for f in os.listdir(directory) ]

        mand = msh.Mesh()
        for f in files:
            if "mand" in f and "warped" in f:
                mand.fondre(msh.Mesh(f))
        mand.write("toto.mesh")
        err1 = os.system( "python ../projects/tools/fastBoundingMesh.py -i toto.mesh -o toto.carved.mesh -r 71")
        err2 = os.system("blender --background --python blender_warp.py -- toto.carved.mesh toto.mesh toto.warped.mesh")
        err3 = os.system("mmgs_O3 toto.warped.mesh -o toto.warped.mesh -hausd 0.001")

    if two:
        cube = msh.Mesh("/home/norgeot/box.mesh")
        mesh = msh.Mesh("toto.warped.mesh")
        print "start"
        dists, _ = nearest_neighbor(cube.verts[:,:3], mesh.verts[:,:3])
        print "end"
        print len(dists), len(mesh.verts), len(cube.verts), dists[0]
        cube.scalars = np.array(dists)
        cube.scaleSol(0.002, 1, absolute=True)
        cube.write("titi.mesh")
        cube.writeSol("titi.sol")
        os.system("mmg3d_O3 titi.mesh -hgrad 1.5")
        os.system("mshdist -noscale titi.o.mesh toto.warped.mesh")
