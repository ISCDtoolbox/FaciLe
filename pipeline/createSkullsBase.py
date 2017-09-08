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

    # 3 - Remesh with two haussdorf distance factors
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


    """
    # 2 - Scale to [0,1]
    MAT = fullMandible.toUnitMatrix()
    np.savetxt("mat_toUnit.txt",MAT)
    fullMandible.applyMatrix(mat=MAT)
    fullMandible.write("mandible.mesh")

    # 4 - Cut the mandible in two
    rightMandible = deepcopy(fullMandible)
    leftMandible = fullMandible
    #Generate the mask
    mask = [1 for i in range(len(leftMandible.tris))]
    mid = np.mean(leftMandible.verts,axis=0)[0]
    print mid
    for i,t in enumerate(leftMandible.tris):
        for v in t:
            x = leftMandible.verts[v][0]
            if x < mid:
                mask[i] = 0
    #Create the left mandible
    leftMandible.tris = np.array([t for i,t in enumerate(leftMandible.tris) if mask[i]==1])
    print len(leftMandible.tris)
    leftMandible.discardUnused()
    leftMAT = leftMandible.toUnitMatrix()
    np.savetxt("1_leftMandibleToUnit.txt", leftMAT)
    leftMandible.applyMatrix(mat=leftMAT)
    leftMandible.write("leftMandible.mesh")
    #And the right one, symetrized
    rightMandible.tris = np.array([t for i,t in enumerate(rightMandible.tris) if mask[i]==0])
    rightMandible.discardUnused()
    rightMAT = rightMandible.toUnitMatrix()
    np.savetxt("1_rightMandibleToUnit.txt", rightMAT)
    rightMandible.applyMatrix(mat=rightMAT)
    rightMandible.verts[:,0] = 1-rightMandible.verts[:,0]
    rightMandible.write("rightMandible.mesh")


    # 5 - Create the shells for left and right mandibles
    #command(exe.boundingMesh + " leftMandible.mesh", displayOutput=True)
    command(exe.shell + " -i leftMandible.mesh -o leftShell.mesh -c", displayOutput=True)
    command(exe.shell + " -i rightMandible.mesh -o rightShell.mesh -c", displayOutput=True)

    sys.exit()


    # 7 - Create a domain for mshdist computation
    #Right mandible
    cube=msh.Mesh(cube=[0,1,0,1,0,1])
    cube.write("rightBox.mesh")
    command( "tetgen -pgANEF rightBox.mesh")
    command( "mmg3d_O3 rightBox.1.mesh -hausd " + str(np.max(mesh.dims)/25) + " -hmax " + str(np.max(mesh.dims)/25))
    command( "mshdist -ncpu 4 -noscale rightBox.1.o.mesh rightMandible.warped.mesh")

    # 9 - Morphing the template_mandibule surface to the computed boxes
    command(morphing + " template_halfMandible_volume.mesh rightBox.1.o.mesh")

    # 10 - Extract the surface from the morphing results
    morphed = msh.Mesh("morphed.mesh")
    morphed.readSol()
    morphed.extractSurfaces#Placeholder
    """
