import os
import sys
import numpy as np
from copy import deepcopy
import argparse

#Parallel
import subprocess as sp
import multiprocessing as mp


sys.path.append(os.path.join(os.path.dirname(__file__),"../projects/tools"))
import msh
import executable_paths as exe

def parse():
    parser = argparse.ArgumentParser(description="Creates mandible and masseter files for the database creation")
    parser.add_argument("-i", "--inputDir", help="input directory", type=str, required=True)
    parser.add_argument("-o", "--outputDir", help="output directory", type=str, required=True)
    return parser.parse_args()
def checkArgs(args):
    if not os.path.exists(args.inputDir):
        print args.input + "is not a valid directory"
        sys.exit()
    if not len([f for f in os.listdir(args.inputDir) if f[0]=="."]) == 0:
        print args.inputDir + " is an empty directory"
        sys.exit()
    if not os.path.exists(args.outputDir):
        print args.outputDir + " does not exist, creating"
        os.system("mkdir " + args.outputDir)
    args.inputDir = os.path.abspath(args.inputDir)
    args.outputDir = os.path.abspath(args.outputDir)
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


def work(in_file):
    """Defines the work unit on an input file"""
    root = '.'.join(in_file.split("/")[-1].split(".")[:-1])

    if not os.path.exists("tmp_"+root):
        os.mkdir("tmp_"+root)
    os.chdir("tmp_"+root)
    os.system("cp /home/norgeot/dev/own/FaciLe/projects/warping/demo/sphere.o1.mesh ./sphere.mesh")

    cmd = " ".join([exe.processSkull, "-i " + in_file, "-t ../../OsTemplate.mesh",">",root+"_OUT.txt"])
    print "Starting the skull processing for " + in_file
    #os.system(cmd)
    print "Skull processing finished for " + in_file

    #clean the working directories

    for ext in [".warped.mesh", ".box.1.o.", "mat","_OUT.txt"]:
        for f in os.listdir("."):
            if ext in f:
                os.rename(f, os.path.join(args.outputDir,f))
    for f in os.listdir("."):
        if ".mesh" in f or ".sol" in f:
            #os.remove(f)
            #print f + " was successfully removed"
            a=2

    return 0

if __name__=="__main__":

    args = parse()
    checkArgs(args)
    files = [os.path.join(args.inputDir,f) for f in os.listdir(args.inputDir) if ".mesh" in f]



    #Set up the parallel task pool to use all available processors
    count = mp.cpu_count()
    pool = mp.Pool(processes=count)
    pool.map(work, files)
