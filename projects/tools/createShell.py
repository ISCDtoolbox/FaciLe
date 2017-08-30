from scipy   import ndimage as nd
from skimage import measure as mea
import numpy as np
import os
import msh
import sys


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
    mesh.tris[:,3]=3



    print "Computing the signed distance"
    cube = msh.Mesh()
    cubeDist = np.max(mesh.dims)/5
    cube.verts, cube.tris = addCube(
        mesh.xmin-cubeDist,
        mesh.xmax+cubeDist,
        mesh.ymin-cubeDist,
        mesh.ymax+cubeDist,
        mesh.zmin-cubeDist,
        mesh.zmax+cubeDist,
        ref=0
    )
    cube.write("box.mesh")
    #Creating the signed distance object
    command( "tetgen -pgANEF box.mesh")
    command( "mmg3d_O3 box.1.mesh -hausd " + str(np.max(mesh.dims)/25) + " -hmax " + str(np.max(mesh.dims)/25))
    command( "mshdist -ncpu 4 -noscale box.1.o.mesh " + sys.argv[1])



    print "Extracting isovalues"
    #Extracting the isovalues meshes
    command("mmg3d_O3 box.1.o.mesh -nr -ls " + str(5) + " -hausd " + str(2) + " -out iso1.mesh")
    iso1 = msh.Mesh("iso1.mesh")
    iso1.tets = np.array([])
    iso1.removeRef(0)
    iso1.discardUnused()
    iso1.write("iso1.mesh")
    #Extracting the isovalues meshes
    command("mmg3d_O3 box.1.o.mesh -nr -ls " + str(1) + " -hausd " + str(2) + " -out iso2.mesh")
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
    final.tris = final.tris[final.tris[:,-1]>0]
    final.tets = final.tets[final.tets[:,-1]!=2]
    final.discardUnused()
    final.write("final.mesh")
    #Last volume remesh
    command("mmg3d_O3 final.mesh -hausd " + str(np.max(final.dims)/1000))



    print "Scaling and translating to a unit box"
    #Scaling to a box between 0 and 1
    shell   = msh.Mesh("final.o.mesh")
    cube    = msh.Mesh(cube=[0.1,0.9,0.1,0.9,0.1,0.9])
    scale  = np.min(cube.dims/shell.dims)
    center = shell.center
    #Applying the scale factors
    print "Scale: " + str(scale)
    shell.scale(scale, center)
    mesh.scale(scale, center)
    #Translating meshes
    print "Translation: " + str(" ".join([ str(x) for x in 0.5-shell.center]) )
    shell.verts[:,:3]+=(np.array([0.5,0.5,0.5])-shell.center)
    mesh.verts[:,:3]+=(np.array([0.5,0.5,0.5])-shell.center)
    shell.computeBBox()
    mesh.computeBBox()
    mesh.tris[:,-1] = 0
    #Writing files
    shell.write("shell.mesh")
    mesh.write("surface.mesh")
