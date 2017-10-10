#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import sys
import numpy as np
from copy import deepcopy
import argparse
from shutil import copyfile

#Parallel
import subprocess as sp
import multiprocessing as mp
from functools import partial

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

def obj2Mesh(file):
    with open(file, "r") as f:
        LINES = f.readlines()
        mesh = msh.Mesh()
        mesh.verts = np.array([ [float(x) for x in l.split()[1:]] for l in LINES if l[0]=="v" ])
        mesh.tris  = np.array([ [int(x)-1 for x in l.split()[1:]]   for l in LINES if l[0]=="f" ])
        mesh.verts = np.insert(mesh.verts,3,0,axis=1)
        mesh.tris  = np.insert(mesh.tris,3,0,axis=1)
        mesh.computeBBox()
        return mesh
        """
        inds = [l.split()[0] for l in LINES]
        nV = len([x for x in inds if x == "v"])
        nT = len([x for x in inds if x == "f"])
        for i,l in enumerate(LINES):
            if (inds[i]=="v"):
                out.write(" ".join([l.split()[i] for i in range(1,4)]) + " 1\n")

        out.write("\nTriangles\n"+str(nT)+"\n")
        for i,l in enumerate(LINES):
            if (inds[i]=="f"):
                out.write(" ".join([str(int(l.split()[i])) for i in range(1,4)]) + " 1\n")

        out.close()
        print "Successfully converted " + inFile
        """
def loadSTL(input, script="cleanSTL.mlx"):
    if os.path.isfile(input):
        name = ".".join(input.split(".")[:-1])
        os.system("LC_ALL=C meshlabserver -i " + input + " -o " + name + ".obj -s " + script + " > /dev/null 2>&1" )
        mesh = obj2Mesh(name+".obj")
        os.remove(name+".obj")
        return mesh
    else:
        print input + " does not exist, exiting..."
        sys.exit()
def group_files(dir):
    files = os.listdir(dir)
    groups = []
    while len(files)>0:
        file = files[0]
        group = [file]
        for other in files:
            if other!=file and other.split(".")[1] == file.split(".")[1]:
                group.append(other)
        groups.append(group)
        for g in group:
            files.remove(g)
    return groups
def newName(f):
    names = [
        ["OsB","bone"],
        ["SkullB","bone"],
        ["MandB","mand"],
        ["BTeeB","btee"],
        ["HTeeB","htee"],
        ["Mass","mass"],
        ["FatB","face"]
    ]
    for n in names:
        if n[0] in f:
            return n[1]

def convertToMesh(inputDir, outputDir, group):
    num = int(group[0].split(".")[1])
    n = str(num).zfill(3)
    des = " ".join(group)

    # 0 - Tout convertir en .mesh
    for f in group:
        try:
            if ".stl" in f:
                mesh = loadSTL(os.path.join(inputDir, f), script="../cleanSTL.mlx")
                mesh.write(os.path.join(outputDir, n + "_" + newName(f) + ".mesh"))
            elif ".mesh" in f:
                copyfile(os.path.join(inputDir, f), os.path.join(outputDir, n + "_" + newName(f) + ".mesh"))
        except:
            print "Error converting " + f

def fondreAndScale(inputDir, outputDir, group):
    des = " ".join(group)

    # face + bone
    if "face" in des and "bone" in des and len(group) == 2:
        boneFile = None

        #Scale and remesh
        for f in group:
            infile  = os.path.join(inputDir, f)
            outfile = os.path.join(outputDir, f)
            mesh = msh.Mesh(infile)
            mesh.verts[:,:3] -= mesh.center
            mesh.verts[:,:3] *= 0.0035
            mesh.verts[:,:3] += [0.5,0.5,0.5]
            mesh.computeBBox()
            mesh.write(outfile)
            if "bone" in f:
                boneFile = os.path.join(outputDir, f[:-5] + ".o.mesh")
            command(exe.mmgs + " " + outfile + " -o " + outfile[:-5] + ".o.mesh -nr -hausd 0.007", displayOutput=True )

        #Align following the skull
        templateFile = "/home/norgeot/dev/own/FaciLe/OsTemplate.mesh"
        command(exe.align + " -i " + boneFile + " " + templateFile + " -d 0.1 -o 0.95", displayOutput=True)
        os.system("mv mat_Super4PCS.txt " + os.path.join(outputDir,f[:3] + "_mat_1_super4pcs.txt"))
        os.system("mv mat_ICP.txt "       + os.path.join(outputDir,f[:3] + "_mat_2_cppICP.txt"))
        for f in group:
            mesh = msh.Mesh(os.path.join(outputDir, f[:-5] + ".o.mesh"))
            mesh.applyMatrix(matFile=os.path.join(outputDir,f[:3] + "_mat_1_super4pcs.txt"))
            mesh.applyMatrix(matFile=os.path.join(outputDir,f[:3] + "_mat_2_cppICP.txt"))
            mesh.write(f[:-5] + ".o.mesh")


if __name__=="__main__":

    args = parse()
    checkArgs(args)

    #Create the directory if it does not exist
    if not os.path.exists(args.outputDir):
        os.mkdir(args.outputDir)
    os.chdir(args.outputDir)

    groups = group_files(args.inputDir)

    # 1 - Move to directory
    """
    count = mp.cpu_count() - 1
    pool = mp.Pool(processes=count)
    func = partial(convertToMesh, args.inputDir, args.outputDir)
    pool.map(func, groups)
    """

    # 2 - Fondre, rescale and align
    newGroups = []
    done      = []
    for f in os.listdir(args.outputDir):
        if f[:3] in done:
            newGroups[done.index(f[:3])].append(f)
        else:
            newGroups.append([f])
            done.append(f[:3])


    count = mp.cpu_count() - 1
    pool = mp.Pool(processes=count)
    func = partial(fondreAndScale, args.outputDir, "/Data/remeshed")
    pool.map(func, newGroups)



    """
    skul = loadSTL("/home/norgeot/LOIC/SkullB.33.stl")
    mand = loadSTL("/home/norgeot/LOIC/MandB.33.stl")
    mand.replaceRef(0,1)
    btee = loadSTL("/home/norgeot/LOIC/BTeeB.33.stl")
    btee.replaceRef(0,2)
    htee = loadSTL("/home/norgeot/LOIC/HTeeB.33.stl")
    htee.replaceRef(0,3)
    mand.fondre(btee)
    mand.fondre(htee)
    skul.fondre(mand)
    skul.write("skull.mesh")
    """

    #Set up the parallel task pool to use all available processors
    #count = mp.cpu_count()
    #pool = mp.Pool(processes=count)
    #pool.map(work, files)
