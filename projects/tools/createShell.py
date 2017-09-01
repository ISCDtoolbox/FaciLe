from scipy   import ndimage as nd
from skimage import measure as mea
import numpy as np
import os
import msh
import sys
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__),"../"))
import executable_paths as exe

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

def parse():
    parser = argparse.ArgumentParser(description="Creates an adequate shell for warping")
    parser.add_argument("-i", "--input", help="input .mesh file", type=str, required=True)
    parser.add_argument("-o", "--output", help="transformed .mesh file", type=str, required=True)
    parser.add_argument("-s", "--scale", help="if specified, scales the input mesh to a [0.2, 0.8] box", action="store_true")
    parser.add_argument("-c", "--carve", help="use a space carve bounding mesh for signed distance computation", action="store_true")
    parser.add_argument("-int", "--interiorIsovalue", help="ratio to the maximum dimension for the interior shell surface distance to the mesh", type=float, default=20)
    parser.add_argument("-ext", "--exteriorIsovalue", help="ratio to the maximum dimension for the exterior shell surface distance to the mesh", type=float, default=10)
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
    if args.interiorIsovalue < args.exteriorIsovalue:
        print "The inner shell must be closer than the outer shell"
        sys.exit()
    if args.interiorIsovalue<5 or args.exteriorIsovalue<5:
        print "The shell must be closer than maxDim/5"
        sys.exit()

if __name__=="__main__":
    args = parse()
    checkArgs(args)

    print "Opening " + args.input
    mesh = msh.Mesh(args.input)

    print "Preparing the domain for signed distance computation"
    if args.scale:
        cube=msh.Mesh(cube=[0,1,0,1,0,1])
        cube.write("box.mesh")
        MAT = mesh.toUnitMatrix()
        np.savetxt(args.input[:-5]+".scalingMat.txt", MAT)
        mesh.applyMatrix(mat=MAT)
        mesh.write(args.input[:-5]+".scaled.mesh")
        if args.carve:
            print "Computing a bounding mesh to use for signed distance computation"
            command(exe.boundingMesh + " -i " + args.input[:-5]+".scaled.mesh" + " -o " + args.input[:-5]+".carved.mesh -r 31")
    else:
        c = np.max(mesh.dims)/8
        cube = msh.Mesh(cube=[mesh.xmin-c,
            mesh.xmax+c,
            mesh.ymin-c,
            mesh.ymax+c,
            mesh.zmin-c,
            mesh.zmax+c])
        cube.write("box.mesh")
        if args.carve:
            print "Computing a bounding mesh to use for signed distance computation"
            command(exe.boundingMesh + " -i " + args.input + " -o " + args.input[:-5]+".carved.mesh -r 31")
    command( "tetgen -pgANEF box.mesh")
    command( "mmg3d_O3 box.1.mesh -hausd " + str(np.max(mesh.dims)/50) + " -hmax " + str(np.max(mesh.dims)/25))

    print "Computing the signed distance"
    import multiprocessing
    ncpus = multiprocessing.cpu_count()
    if args.carve:
        command( "mshdist -ncpu " + str(ncpus) + " -noscale box.1.o.mesh " + args.input[:-5]+".carved.mesh")
    else:
        if args.scale:
            command( "mshdist -ncpu " + str(ncpus) + " -noscale box.1.o.mesh " + args.input[:-5]+".scaled.mesh")
        else:
            command( "mshdist -ncpu " + str(ncpus) + " -noscale box.1.o.mesh " + args.input)



    print "Extracting isovalues"
    #Extracting the isovalues meshes
    command("mmg3d_O3 box.1.o.mesh -nr -ls " + str(np.max(mesh.dims)/args.exteriorIsovalue) + " -hausd " + str(np.max(mesh.dims)/50) + " -out iso1.mesh")
    iso1 = msh.Mesh("iso1.mesh")
    iso1.tets = np.array([])
    iso1.removeRef(0)
    iso1.discardUnused()
    iso1.write("iso1.mesh")
    #Extracting the isovalues meshes
    command("mmg3d_O3 box.1.o.mesh -nr -ls " + str(np.max(mesh.dims)/args.interiorIsovalue) + " -hausd " + str(np.max(mesh.dims)/50) + " -out iso2.mesh")
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
    final.write(args.output)
    #Last volume remesh
    command("mmg3d_O3 " + args.output +  " -hausd " + str(np.max(final.dims)/1000))
