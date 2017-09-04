from scipy   import ndimage as nd
from skimage import measure as mea
import numpy as np
import os
import msh
import sys
import argparse
sys.path.append(os.path.join(os.path.dirname(__file__),"../"))
import executable_paths as exe

def parse():
    parser = argparse.ArgumentParser(description="Creates an adequate shell for warping")
    parser.add_argument("-i", "--input", help="input .mesh file", type=str, required=True)
    parser.add_argument("-o", "--output", help="transformed .mesh file", type=str, required=True)
    parser.add_argument("-c", "--center", help="center of the icosphere", type=float, default=[0.,0.,0.], metavar=('center x', 'center y', 'center z'), nargs=3)
    parser.add_argument("-r", "--radius", help="radius of the icosphere", type=float, default=0.1)
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
    lines = os.popen("tetgen -pdNEF " + args.input).read()
    if "intersects" in lines:
        print args.input + " has intersecting facets"
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
    """
    Suppose we are between 0 and 1 for a skull and take a wrapped mesh
    """

    args = parse()
    checkArgs(args)

    #command("mmgs_O3 " + sys.argv[1] + " -nr -o out.mesh")
    #mesh = msh.Mesh("out.mesh")
    mesh = msh.Mesh(args.input)
    ico  = msh.Mesh(ico=[args.center,args.radius])
    ico.tris[:,-1]=10
    mesh.fondre(ico)
    mesh.write("out.mesh")
    command("tetgen -pgANEF out.mesh")
    command("mmg3d_O3 out.1.mesh -nosurf -o " + args.output)
    command("rm out.mesh out.1.mesh")
