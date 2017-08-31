import msh
import sys

if __name__=="__main__":
    mesh = msh.Mesh(sys.argv[1])
    mesh.applyMatrix(matFile=sys.argv[2])
    mesh.write("out.mesh")
