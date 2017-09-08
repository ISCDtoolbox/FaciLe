import msh
import sys
import argparse
import os
import numpy as np

def parse():
    parser = argparse.ArgumentParser(description="Computes the distances between two .mesh files")
    parser.add_argument("-i1", "--input1", help="input .mesh file", type=str, required=True)
    parser.add_argument("-i2", "--input2", help="input .mesh file", type=str, required=True)
    parser.add_argument("-o", "--output", help="output.sol file", type=str, required=True)
    return parser.parse_args()

def checkArgs(args):
    if not os.path.isfile(args.input1):
        print args.input1 + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.input1)[1] == ".mesh":
        print args.input1 + " is not a .mesh file"
        sys.exit()
    if not os.path.isfile(args.input2):
        print args.input2 + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.input2)[1] == ".mesh":
        print args.input2 + " is not a .mesh file"
        sys.exit()
    if not os.path.splitext(args.output)[1] == ".mesh" and not os.path.splitext(args.output)[1] == ".sol":
        print "Output file must be in the .mesh format"
        sys.exit()

if __name__=="__main__":
    args = parse()
    checkArgs(args)
    mesh1 = msh.Mesh(args.input1)
    mesh2 = msh.Mesh(args.input2)
    if len(mesh1.verts) != len(mesh2.verts):
        print "Different number of verts for the two meshes"
        sys.exit()
    mesh1.scalars = np.linalg.norm(mesh1.verts[:,:3] - mesh2.verts[:,:3], axis=1)
    if os.path.splitext(args.output)[1] == ".mesh":
        mesh1.write(args.output)
        mesh1.writeSol(args.output[:-5] + ".sol")
    else:
        mesh1.writeSol(args.output)
