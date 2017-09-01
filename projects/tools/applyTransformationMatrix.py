import msh
import sys
import argparse
import os

def parse():
    parser = argparse.ArgumentParser(description="Apply a transformation matrix to a .mesh file")
    parser.add_argument("-i", "--input", help="input .mesh file", type=str, required=True)
    parser.add_argument("-m", "--matrix", help="input matrix file", type=str, required=True)
    parser.add_argument("-o", "--output", help="transformed .mesh file", type=str, required=True)
    return parser.parse_args()

def checkArgs(args):
    if not os.path.isfile(args.input):
        print args.input + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.input)[1] == ".mesh":
        print args.input + " is not a .mesh file"
        sys.exit()
    if not os.path.isfile(args.matrix):
        print args.matrix + " is not a valid file"
        sys.exit()
    if not os.path.splitext(args.output)[1] == ".mesh":
        print "Output file must be in the .mesh format"
        sys.exit()

if __name__=="__main__":
    args = parse()
    checkArgs(args)
    mesh = msh.Mesh(args.input)
    mesh.applyMatrix(matFile=args.matrix)
    mesh.write(args.output)
