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
    args.template = os.path.abspath(args.template)
    args.input    = os.path.abspath(args.input)
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
    print root

    #1 - Load the mesh
    skull = msh.Mesh(args.input)

    # 2 - Scale by a fixed factor and align to the template
    skull.verts[:,:3] -= skull.center
    skull.verts[:,:3] *= 0.0035
    skull.verts[:,:3] += [0.5,0.5,0.5]
    skull.computeBBox()
    skull.write(root + ".scaled.mesh")

    # 3 - Remesh with two hausdorff distance factors
    command(exe.mmgs + " " + root + ".scaled.mesh -o " + root + ".o1.mesh -nr -hausd 0.007" )
    command(exe.mmgs + " " + root + ".scaled.mesh -o " + root + ".o2.mesh -nr -hausd 0.0007" )

    # 4 - Align to the template skull
    command(exe.align + " -i " + root + ".o1.mesh "+args.template+" -d 0.1 -o 0.95", displayOutput=True)
    os.system("mv mat_Super4PCS.txt " + root + "_mat_1_super4pcs.txt")
    os.system("mv mat_ICP.txt " + root + "_mat_2_cppICP.txt")
    command(exe.pythonICP + " -s " + root + ".o1.mesh -t " + args.template + " -m "+root+"_mat_3_pyICP.txt")

    # 5 - Apply the transformation matrices
    remeshed = msh.Mesh(root + ".o1.mesh")
    remeshed.applyMatrix(matFile=root + "_mat_1_super4pcs.txt")
    remeshed.applyMatrix(matFile=root + "_mat_2_cppICP.txt")
    remeshed.applyMatrix(matFile=root + "_mat_3_pyICP.txt")
    remeshed.write(root + ".o1.mesh")

    # 6 - Create the warping shell with carving
    #command(exe.shell + " -i " + root + ".o1.mesh -o shell.mesh -c")

    # 7 - Warp
    #command(exe.warping + " " + root + ".o1.mesh -t shell.mesh -p -load 50", displayOutput=True)
    command(exe.warping + " " + root + ".o1.mesh -p -load 70", displayOutput=True)

    # 8 - Extract the interior surface
    warped = msh.Mesh("sphere.d.mesh")
    ext_ref = 2
    warped.tris = warped.tris[warped.tris[:,-1] != ext_ref]
    warped.tets = np.array([])
    warped.discardUnused()
    warped.write(root + ".warped.mesh")

    # 9 - ICP alignement on the warped meshes
    command(exe.pythonICP + " -s " + root + ".warped.mesh -t " + args.template + " -m "+root+"_mat_4_pyICP.txt -mIts 300 -tol 0.00001 -mPts 10000")
    warped.applyMatrix(matFile=root + "_mat_4_pyICP.txt")
    warped.write(root + ".warped.mesh")

    # 10 - Signed distance
    cube=msh.Mesh(cube=[0,1,0,1,0,1])
    cube.write(root + ".box.mesh")
    command( exe.tetgen + " -pgANEF " + root + ".box.mesh")
    command( exe.mmg3d + " "+root+".box.1.mesh -hausd " + str(np.max(remeshed.dims)/25) + " -hmax " + str(np.max(remeshed.dims)/25))
    command( exe.mshdist + " -ncpu 4 -noscale "+root+".box.1.o.mesh " + root + ".warped.mesh")


    # 11 - Morph the reference onto the skull and extract the surface displacement
    """
    command(exe.morphing + "box.1.o.mesh template_skull.mesh")
    #To do with chiara's code later, here at least it works
    template = msh.Mesh("template_skull.mesh")
    morphed  = msh.Mesh("morphed.mesh")
    dists = [ np.linalg.norm(v1-v2) for v1,v2 in zip(template.verts[:,:3], morphed.verts[:,:3]) ]
    """
    # 10 - Generate "la masque"
