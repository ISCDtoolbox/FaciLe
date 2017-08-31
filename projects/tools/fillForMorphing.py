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
    """
    Suppose we are between 0 and 1 for a skull and take a wrapped mesh
    """
    #command("mmgs_O3 " + sys.argv[1] + " -nr -o out.mesh")
    #mesh = msh.Mesh("out.mesh")
    mesh = msh.Mesh(sys.argv[1])
    ico  = msh.Mesh(ico=[[0.5,0.5,0.6],np.max(mesh.dims)/8])
    ico.tris[:,-1]=10
    mesh.fondre(ico)
    mesh.write("out.mesh")
    command("tetgen -pgANEF out.mesh")
    command("mmg3d_O3 out.1.mesh -nosurf -o out.mesh")
    command("rm out.1.mesh out.sol")
