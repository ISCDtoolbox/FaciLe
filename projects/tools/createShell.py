from scipy   import ndimage as nd
from skimage import measure as mea
import numpy as np
import os
import msh
import sys

def command(cmd, displayOutput=False):
    err = 1
    print "Running the command '" + cmd + "'"

    if displayOutput:
        err = os.system(cmd)
    else:
        err = os.system(cmd + " > tmp_out.txt 2>tmp_err.txt")

    if err:
        print "An error happened while executing:\n"+cmd+"\nLook in tmp_out.txt or tmp_err.txt for info\nExiting..."
        sys.exit()
    else:
        os.system("rm tmp_out.txt tmp_err.txt >/dev/null 2>&1")

if __name__=="__main__":


    print "Opening " + sys.argv[1]
    mesh = msh.Mesh(sys.argv[1])

    print "Computing the signed distance"
    """
    cubeDist = np.max(mesh.dims)/5
    cube = msh.Mesh(cube=[mesh.xmin-cubeDist,
        mesh.xmax+cubeDist,
        mesh.ymin-cubeDist,
        mesh.ymax+cubeDist,
        mesh.zmin-cubeDist,
        mesh.zmax+cubeDist])
    """
    cube=msh.Mesh(cube=[0,1,0,1,0,1])
    cube.write("box.mesh")
    #Creating the signed distance object
    command( "tetgen -pgANEF box.mesh")
    command( "mmg3d_O3 box.1.mesh -hausd " + str(np.max(mesh.dims)/25) + " -hmax " + str(np.max(mesh.dims)/25))
    command( "mshdist -ncpu 4 -noscale box.1.o.mesh " + sys.argv[1])



    print "Extracting isovalues"
    #Extracting the isovalues meshes
    command("mmg3d_O3 box.1.o.mesh -nr -ls " + str(np.max(mesh.dims)/10) + " -hausd " + str(np.max(mesh.dims)/50) + " -out iso1.mesh")
    iso1 = msh.Mesh("iso1.mesh")
    iso1.tets = np.array([])
    iso1.removeRef(0)
    iso1.discardUnused()
    iso1.write("iso1.mesh")
    #Extracting the isovalues meshes
    command("mmg3d_O3 box.1.o.mesh -nr -ls " + str(np.max(mesh.dims)/20) + " -hausd " + str(np.max(mesh.dims)/50) + " -out iso2.mesh")
    iso2 = msh.Mesh("iso2.mesh")
    iso2.tets = np.array([])
    iso2.removeRef(0)
    iso2.discardUnused()
    iso2.write("iso2.mesh")
    #Merging meshes
    iso1.replaceRef(10,1)
    iso2.replaceRef(10,2)
    iso1.fondre(iso2)
    iso1.write("iso3.mesh")



    print "Creating a volume mesh for the shell"
    #Creating tetrahedra
    command("tetgen -pgaYANEFq1.02 iso3.mesh")
    final = msh.Mesh("iso3.1.mesh")

    #Pour un triangle de l'exterieur (ref 1), si un tetra partage un point,
    #alors on garde sa reference
    ext_point_ind = final.tris[final.tris[:,-1]==1][0][0]
    ext_ref=None
    for t in final.tets:
        if ext_point_ind in t:
            ext_ref = t[-1]

    final.tris = final.tris[final.tris[:,-1]>0]
    final.tets = final.tets[final.tets[:,-1]==ext_ref]
    final.discardUnused()
    final.write("final.mesh")
    #Last volume remesh
    command("mmg3d_O3 final.mesh -hausd " + str(np.max(final.dims)/1000))



    print "Scaling and translating to a unit box"
    """
    #Scaling to a box between 0 and 1
    shell   = msh.Mesh("final.o.mesh")
    scale  = 0.8/np.max(shell.dims)
    #Translate to origin
    M1 = np.eye(4)
    M1[:3,3] = -shell.center
    #Scale
    M2 = np.eye(4)*scale
    M2[3,3]=1
    #Translate to 0.5
    M3 = np.eye(4)
    M3[:3,3] = np.array([0.5,0.5,0.5])
    MAT = np.dot(np.dot(M3, M2),M1)
    np.savetxt("matToUnit.txt",MAT)
    """
