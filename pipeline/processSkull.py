import os
import sys
import numpy as np
from copy import deepcopy
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__),"../projects/tools"))
import msh
import executable_paths as exe

def parse():
    parser = argparse.ArgumentParser(description="Creates mandible and masseter files for the database creation")
    parser.add_argument("-i", "--input", help="input .mesh object", type=str, required=True)
    parser.add_argument("-t", "--template", help="template mandible in full size", type=str, required=True)
    parser.add_argument("-v", "--verbose", help="talk to me", action="store_true")
    return parser.parse_args()

def checkArgs(args):
    if not os.path.isfile(args.input):
        print args.input + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.input)[1] == ".mesh":
        print args.input + " is not a .mesh file"
        sys.exit()
    if not os.path.isfile(args.template):
        print args.template + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.template)[1] == ".mesh":
        print args.template + " is not a .mesh file"
        sys.exit()

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

    args = parse()
    checkArgs(args)
    root = ".".join(args.input.split("/")[-1].split(".")[:-1])
    """

    #1 - Load the mesh
    skull = msh.Mesh(args.input)

    # 2 - Scale by a fixed factor
    skull.verts[:,:3] -= skull.center
    skull.verts[:,:3] *= 0.0035
    skull.verts[:,:3] += [0.5,0.5,0.5]
    skull.computeBBox()
    skull.write(root + ".scaled.mesh")


    # 3 - Remesh with two haussdorf distance factors
    command(exe.mmgs + " " + root + ".scaled.mesh -o " + root + ".scaled.o1.mesh -nr -hausd " + str(np.max(skull.dims)/100.) )
    command(exe.mmgs + " " + root + ".scaled.mesh -o " + root + ".scaled.o2.mesh -nr -hausd " + str(np.max(skull.dims)/1000.) )

    # 4 - Align to the template skull
    command(exe.align + " -i " + root + ".scaled.o2.mesh " + args.template + " -d 0.1 -o 0.85", displayOutput=True)
    #command(exe.pythonICP + " -s " + root + ".scaled.o1.mesh -t " + args.template + " -m mat_PythonICP.txt")

    # 5 - Applying the transformation matrices
    skull1 = msh.Mesh(root + ".scaled.o1.mesh")
    skull2 = msh.Mesh(root + ".scaled.o2.mesh")
    objs = [skull, skull1, skull2]
    names = [".scaled.mesh", ".scaled.o1.mesh", ".scaled.o2.mesh"]
    #objs = [skull1, skull2]
    #names = [".scaled.o1.mesh", ".scaled.o2.mesh"]
    for obj,ext in zip(objs, names):
        obj.applyMatrix(matFile="mat_Super4PCS.txt")
        obj.applyMatrix(matFile="mat_ICP.txt")
        #obj.applyMatrix(matFile="mat_PythonICP.txt")
        obj.write(root + ext)


    # 6 - Create the warping shell with carving
    command(exe.shell + " -i " + root + ".scaled.o2.mesh -o " + root + ".shell.mesh -c", displayOutput=True)
    """


    # 7 - Warp
    command(exe.warping + " -t " + root + ".shell.o.mesh -s " + root + ".scaled.o2.mesh", displayOutput=True)
    sys.exit()

    # 8 - Signed distance
    cube=msh.Mesh(cube=[0,1,0,1,0,1])
    cube.write("box.mesh")
    command( "tetgen -pgANEF box.mesh")
    command( "mmg3d_O3 box.1.mesh -hausd " + str(np.max(skull.dims)/25) + " -hmax " + str(np.max(skull.dims)/25))
    command( "mshdist -ncpu 4 -noscale box.1.o.mesh skull.warped.mesh")

    # 8 - Morph the reference onto the skull and extract the surface displacement
    command(exe.morphing + " template_skull.mesh box.1.o.mesh")
    #To do with chiara's code later, here at least it works
    template = msh.Mesh("template_skull.mesh")
    morphed  = msh.Mesh("morphed.mesh")
    dists = [ np.linalg.norm(v1-v2) for v1,v2 in zip(template.verts[:,:3], morphed.verts[:,:3]) ]

    # 9 - Generate "la masque"
    #???
