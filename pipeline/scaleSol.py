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
def adapt_box_to(f, maxNb=30000):
    cube = msh.Mesh("/home/norgeot/box.mesh")
    mesh = msh.Mesh(f)
    step = 1 if len(mesh.verts)<maxNb else int(len(mesh.verts)/maxNumPoints)+1
    dists, _ = nearest_neighbor(cube.verts[:,:3], mesh.verts[::step,:3])
    cube.scalars = np.array(dists)
    cube.scaleSol(0.002, 1, absolute=True)
    cube.write("cube.mesh")
    cube.writeSol("cube.sol")
    err = os.system("mmg3d_O3 cube.mesh -o cube.o.mesh -hgrad 1.5")
    if err:
        raise IOError("mmg3d failure")

if __name__=="__main__":
    f = "/Data/Facile2/muscles/014_mand.R.warped.mesh"
    adapt_box_to(f)
    os.system("mshdist -noscale cube.o.mesh " + f)

    """
    if one:
        directory = "/Data/Facile2/remeshed/"
        files = [ os.path.join(directory, f) for f in os.listdir(directory) ]

        mesh = msh.Mesh()
        for f in files:
            if "face" in f and "final.o" in f:
                mesh.fondre(msh.Mesh(f))
        mesh.write("toto.mesh")
        err1 = os.system("python ../projects/tools/fastBoundingMesh.py -i toto.mesh -o toto.carved.mesh -r 71")
        err2 = os.system("blender --background --python blender_warp.py -- toto.carved.mesh toto.mesh toto.warped.mesh")
        err3 = os.system("mmgs_O3 toto.warped.mesh -o toto.warped.mesh -hausd 0.001")
    """
