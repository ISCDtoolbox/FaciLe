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

    #1 - Load the mesh
    skull = msh.Mesh(args.input)

    # 2 - Scale by a fixed factor and align to the template
    skull.verts[:,:3] *= 0.0035
    skull.verts[:,:3] += [0.5,0.5,0.5] - skull.center
    skull.computeBBox()
    skull.write("skull.mesh")

    # 3 - Remesh with two hausdorff distance factors
    command(exe.mmgs + " skull.mesh -o skull.o1.mesh -nr -hausd " + str(np.max(skull.dims)/100.) )
    command(exe.mmgs + " skull.mesh -o skull.o2.mesh -nr -hausd " + str(np.max(skull.dims)/1000.) )

    # 4 - Align to the template skull
    command(exe.align + " -i skull.o1.mesh -d 50 -o 0.95", displayOutput=True)
    command(exe.pythonICP + " -s skull.o1.mesh -t " + args.template + " -m mat_PythonICP.txt")
    #fullMandible.applyMatrix(matFile="mat_Super4PCS.txt")
    #fullMandible.applyMatrix(matFile="mat_PythonICP.txt")

    # 5 -Create the warping shell with carving
    command(exe.shell + " - i " + args.input[:-5] + ".o1.mesh -o shell.mesh -c")

    # 6 - Warp
    command(exe.warping + " shell.mesh skull.o2.mesh")

    # 7 - Signed distance
    cube=msh.Mesh(cube=[0,1,0,1,0,1])
    cube.write("box.mesh")
    command( "tetgen -pgANEF box.mesh")
    command( "mmg3d_O3 box.1.mesh -hausd " + str(np.max(skull.dims)/25) + " -hmax " + str(np.max(skull.dims)/25))
    command( "mshdist -ncpu 4 -noscale box.1.o.mesh skull.warped.mesh")
    
    # 8 - Create the reference skull 
    
    # 9 - Morph the reference onto the skull and extract the surface displacement
    command(exe.morphing + "box.1.o.mesh template_skull.mesh")
    #To do with chiara's code later, here at least it works
    template = msh.Mesh("template_skull.mesh")
    morphed  = msh.Mesh("morphed.mesh")
    dists = [ np.linalg.norm(v1-v2) for v1,v2 in zip(template.verts[:,:3], morphed.verts[:,:3]) ]

    # 10 - Generate "la masque"
    #???
