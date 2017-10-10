#!/usr/bin/python
# -*- coding: utf-8 -*

import os
import glob
import sys
import numpy as np
import getpass
from ftplib import FTP
import shutil
import subprocess as sp
import multiprocessing as mp
sys.path.append(os.path.join(os.path.dirname(__file__),"../projects/tools"))
import msh
import executable_paths as exe

#Functions to group files by name
def condition1(ref, tmp):
    #Return 1 if the number splitted between dots are the same
    return ref != tmp and ref.split(".")[1] == tmp.split(".")[1]
def condition2(ref, tmp):
    #Return 1 if the first 3 letters are the same
    return ref != tmp and ref[:3] == tmp[:3]
def group(files, conditionFunc):
    groups = []
    while len(files)>0:
        group = [files[0]]
        for other in files:
            if conditionFunc(files[0], other):
                group.append(other)
        groups.append(group)
        for g in group:
            files.remove(g)
    return groups

#Functions to copy files and change their names
def readCSV(filename):
    good = []
    with open(filename) as f:
        LINES = [l.strip().split(",") for l in f.readlines()]
        for i,l in enumerate(LINES):
            #A fermer = on ne traite pas
            if "fermer" in l[4] or "NON" in l[4]:
                good.append(0)
            #inférieur à 24
            elif i<25:
                good.append(0)
            else:
                good.append(1)
        return good
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

#Functions to load and process .obj and .stl
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
def stlToMesh(f):
    try:
        tmp = f[:-4] + ".obj"
        os.system("LC_ALL=C meshlabserver -i " + f + " -o " + tmp + " -s " + mshScript + " > /dev/null 2>&1" )
        mesh = obj2Mesh(tmp)
        mesh.write(f[:-4] + ".mesh")
        os.remove(tmp)
        os.remove(f)
        return 1
    except:
        return -1
def cutMeshInHalf(mesh):
    half1 = msh.Mesh()
    half2 = msh.Mesh()
    half1.verts = mesh.verts
    half2.verts = mesh.verts

    mid = (np.max(mesh.verts[:,0]) + np.min(mesh.verts[:,0]))/2.
    vertsMask = [v[0]<mid for v in mesh.verts]

    mask = [1 for i in range(len(mesh.tris))]
    for i,t in enumerate(mesh.tris):
        for v in t[:3]:
            mask[i] = vertsMask[v]

    half1.tris  = np.array([t for i,t in enumerate(mesh.tris) if mask[i]])
    half2.tris  = np.array([t for i,t in enumerate(mesh.tris) if not mask[i]])
    del mesh

    half1.discardUnused()
    half2.discardUnused()
    half2.verts[:,0] = 1-half2.verts[:,0]
    half2.computeBBox()
    half1.applyMatrix(half1.toUnitMatrix())
    half2.applyMatrix(half2.toUnitMatrix())
    return half1, half2

#Functions to run in parallel on all files
dryRun = False
def run(func, liste, parallel=True):
    if not dryRun:
        if parallel:
            pool = mp.Pool(processes= min(len(liste), mp.cpu_count()-1 ) )
            return pool.map(func, liste )
        else:
            res = []
            for l in liste:
                res.append(func(l))
            return res
    else:
        return 0
def copy(f):
    ftp = FTP(IPadress, ftpUsr, ftpPwd)
    ftp.cwd(ftpDir)
    num = f.split(".")[1].zfill(3)
    localFile = os.path.join(outFolder, num + "_" + newName(f) + "." + f.split(".")[-1])
    if not num + "_" + newName(f) + ".mesh" in "".join(os.listdir(outFolder)):
        if not os.path.isfile(localFile):
            print "Copying from ftp to " + localFile
            with open(localFile, 'wb') as ff:
                ftp.retrbinary('RETR %s' % f, ff.write)
def convertToMesh(f):
    if "stl" in f:
        print "Converting " + f + " to .mesh"
        if not stlToMesh(f):
            print "-> Error converting " + f + " to .mesh"
def cleanMesh(f):
    print "Analysing " + f
    mesh = msh.Mesh()
    mesh.get_infos(f)
    nV = mesh.numItems[0]
    nT = mesh.numItems[1]
    if nV>nT:
        print "Cleaning " + f
        mesh = msh.Mesh(f)
        mesh.discardDuplicateVertices()
        mesh.discardUnused()
        mesh.write(f)
        print "Done cleaning " + f
def process(g):
    newBone = os.path.join(remeshedFolder, g[0][:3] + "_bone.mesh")
    newFace = os.path.join(remeshedFolder, g[0][:3] + "_face.mesh")
    newMass = os.path.join(remeshedFolder, g[0][:3] + "_mass.mesh")
    print "Processing group " + str(g)
    faceFile = [f for f in g if "face" in f][0]
    face = msh.Mesh(os.path.join(outFolder, faceFile))
    boneFiles = [f for f in g if f!=faceFile and "mass" not in f]
    bone   = None
    center = [0,0,0]
    bone = msh.Mesh(os.path.join(outFolder,boneFiles[0]))
    if len(boneFiles)>1:
        for f in boneFiles[1:]:
            bone.fondre(msh.Mesh(os.path.join(outFolder,f)))
    center   = bone.center
    scale    = 0.0035
    meshes = [face, bone]
    mass = None
    if "mass" in "".join(g):
        mass = msh.Mesh(os.path.join(outFolder, [f for f in g if "mass" in f][0]))
        meshes.append(mass)
    for mesh in meshes:
        mesh.verts[:,:3] -= center
        mesh.verts[:,:3] *= 0.0035
        mesh.verts[:,:3] += [0.5,0.5,0.5]
        mesh.computeBBox()
    bone.write(newBone)
    face.write(newFace)
    if mass is not None:
        mass.write(newMass)
def remesh(f, hausd=0.005):
    if not os.path.isfile(f[:-5] + ".o.mesh"):
        print "Remeshing " + f
        err = os.system(exe.mmgs + " " + f + " -nr -nreg -hausd " + str(hausd) + " -o " + f[:-5] + ".o.mesh > " + f[:-5] + ".remesh.txt 2>&1")
        if err:
            print "Error while remeshing " + f + ", look in " + f[:-5] + ".remesh.txt"
            return 1
        else:
            os.remove(f[:-5] + ".remesh.txt")
            print "Successfully remeshed " + f
            return 0
    return 1
def align(g):
    boneFile = [f for f in g if "bone" in f][0]
    faceFile = [f for f in g if "face" in f][0]
    massFile = [f for f in g if "mass" in f][0] if "mass" in "".join(g) else ""
    num = boneFile.split("/")[-1][:3]
    name = "".join(f.split("/")[-1].split(".")[:-1])
    try:
        os.makedirs(name)
    except:
        pass
    os.chdir(name)

    print "Aligning " + boneFile
    err = os.system(exe.align + " -i " + boneFile + " " + boneTemplate + " -d 0.1 -o 0.95  > "+num+".txt")#/dev/null 2>&1")
    if err:
        print "-- Error aligning " + boneFile
    else:
        bone = msh.Mesh(boneFile)
        bone.applyMatrix(matFile = "mat_Super4PCS.txt")
        bone.applyMatrix(matFile = "mat_ICP.txt")
        bone.write(boneFile[:-5]+".aligned.mesh")
        print num + " bone written"

        err = os.system(exe.pythonICP + " -s " + boneFile[:-5]+".aligned.mesh" + " -t " + boneTemplate + " -m mat_pyICP.txt >> " + num + ".txt")
        if err:
            print "-- Error with ICP for " + boneFile
        else:
            bone.applyMatrix(matFile = "mat_pyICP.txt")
            bone.write(boneFile[:-5]+".final.mesh")
            print num + " final bone written"

        face = msh.Mesh(faceFile)
        face.applyMatrix(matFile = "mat_Super4PCS.txt")
        face.applyMatrix(matFile = "mat_ICP.txt")
        face.applyMatrix(matFile = "mat_pyICP.txt")
        face.write(faceFile[:-5]+".final.mesh")
        print num + " final face written"

        if ".mesh" in massFile:
            mass = msh.Mesh(massFile)
            mass.applyMatrix(matFile = "mat_Super4PCS.txt")
            mass.applyMatrix(matFile = "mat_ICP.txt")
            mass.applyMatrix(matFile = "mat_pyICP.txt")
            mass.write(massFile[:-5]+".final.mesh")
            print num + " final mass written"

    os.chdir("..")
    return 0
def warp(f):
    num = f.split("/")[-1][:3]
    name = "".join(f.split("/")[-1].split(".")[:-1])
    if not os.path.exists(name):
        os.makedirs(name)
    os.chdir(name)
    os.system("cp " + templateSphere + " ./sphere.mesh")
    print "Warping " + f
    err = os.system( exe.warping + " " + f + " -p -load 70 > warping.txt" )
    if err:
        print "-- Error while wrapping " + f + ", trying to remove elements"
        try:
            mesh = msh.Mesh(f)
            mesh.tris = np.array([ t for t in mesh.tris if np.max(t)<=len(mesh.verts) ])
            mesh.discardUnused()
            mesh.write(f)
            err = os.system( exe.warping + " " + f + " -p -load 70 > warping.txt" )
            if err:
                os.chdir("..")
                return 1
        except:
            print "Problem not solved"
            os.chdir("..")
            return 1

    warped = msh.Mesh("sphere.d.mesh")
    ext_ref = 2
    warped.tris = warped.tris[warped.tris[:,-1] != ext_ref]
    warped.tets = np.array([])
    warped.discardUnused()
    warped.write(".".join(f.split(".")[:-1]) + ".warped.mesh")
    print "Successfully warped " + f

    os.chdir("..")
    return 0
def signedDistance(f):
    print "Computing the signed distance for " + f
    num = f.split("/")[-1][:3]
    name = "".join(f.split("/")[-1].split(".")[:-1])
    if not os.path.exists(name):
        os.makedirs(name)
    os.chdir(name)

    cube=msh.Mesh(cube=[0,1,0,1,0,1])
    cube.write("box.mesh")

    err = os.system( exe.tetgen + " -pgANEF box.mesh > error.txt")
    if err:
        print "Error with tetgen on " + f
        os.chdir("..")
        return 1
    err = os.system( exe.mmg3d + " box.1.mesh -hausd 0.04 -hmax 0.04 >> error.txt" )
    if err:
        print "Error with mmg3d on " + f
        os.chdir("..")
        return 1
    err = os.system( exe.mshdist + " -ncpu 1 -noscale box.1.o.mesh " + f + " >> error.txt")
    if err:
        print "Error with mshdist on " + f
        os.chdir("..")
        return 1

    os.system("mv box.1.o.mesh " + ".".join(f.split(".")[:-2]) + ".signed.mesh")
    os.system("mv box.1.o.sol  " + ".".join(f.split(".")[:-2]) + ".signed.sol")

    print "Successfully computed the signed distance for " + f

    os.chdir("..")
    shutil.rmtree(name)
    return 0
def cut(f):
    print "Cutting " + f
    mesh=msh.Mesh(os.path.join(outFolder, f))
    MAT = mesh.toUnitMatrix()
    mesh.applyMatrix(MAT)
    newF = os.path.join(musclesFolder, f[:-5] + ".scaled.mesh")
    mesh.write(newF)
    err = remesh(newF, hausd=0.001)
    mesh = msh.Mesh(os.path.join(musclesFolder, f[:-5] + ".scaled.o.mesh"))
    half1, half2 = cutMeshInHalf(mesh)
    half1.write(os.path.join(musclesFolder, f[:-5]+".R.mesh"))
    half2.write(os.path.join(musclesFolder, f[:-5]+".L.mesh"))
    os.remove(os.path.join(musclesFolder, f[:-5] + ".scaled.mesh"))
    os.remove(os.path.join(musclesFolder, f[:-5] + ".scaled.o.mesh"))
    os.remove(os.path.join(musclesFolder, f[:-5] + ".scaled.o.sol"))
    print "Successfully cut " + f
def alignToTemplate(f):
    num = f.split("/")[-1][:3]
    try:
        os.makedirs(num+"align")
    except:
        pass
    os.chdir(num+"align")

    print "Aligning " + f
    err=0
    """
    if "mass" in f:
        err = os.system(exe.align + " -i " + f + " " + massTemplate + " -d 0.1 -o 0.95  > "+num+"_mass.txt")#/dev/null 2>&1")
    elif "mand" in f:
        err = os.system(exe.align + " -i " + f + " " + mandTemplate + " -d 0.1 -o 0.95  > "+num+"_mand.txt")#/dev/null 2>&1")

    if err:
        print "-- Error aligning " + boneFile
    """
    if 0:
        pass
    else:
        mesh = msh.Mesh(f)
        #mesh.applyMatrix(matFile = "mat_Super4PCS.txt")
        #mesh.applyMatrix(matFile = "mat_ICP.txt")
        #mesh.write(f[:-5]+".aligned.mesh")
        todo = f#f[:-5]+".aligned.mesh"
        if "mass" in f:
            err = os.system(exe.pythonICP + " -s " + todo + " -t " + massTemplate + " -m mat_pyICP.txt >> " + num + "_mass.txt")
        elif "mand" in f:
            err = os.system(exe.pythonICP + " -s " + todo + " -t " + mandTemplate + " -m mat_pyICP.txt >> " + num + "_mand.txt")
        if err:
            print "-- Error with ICP for " + f
        else:
            mesh.applyMatrix(matFile = "mat_pyICP.txt")
            mesh.write(f)
            print "Successfully aligned " + f

    os.chdir("..")
    return 0
def warpWithBlender(f,resolution=41):
    print "Warping "+f+" via blender"
    err1 = os.system( exe.boundingMesh + " -i " + f + " -o " + f[:-5]+".carved.mesh -r " + str(resolution) + " > /dev/null 2>&1")
    err2 = os.system("blender --background --python blender_warp.py -- " + f[:-5]+".carved.mesh " + f + " " + f[:-5] + ".warped.mesh > /dev/null 2>&1")
    err3 = os.system(exe.mmgs + " " + f[:-5] + ".warped.mesh -o " + f[:-5] + ".warped.mesh -hausd 0.002 > /dev/null 2>&1")
    if (err1+err2+err3):
        print "Did not manage to warp "+f+" with blender "
    else:
        os.remove(f[:-5]+".carved.mesh")
        os.remove(f[:-5]+".warped.sol")
def generateMask(f):
    num = f[0].split("/")[-1][:3]
    try:
        os.makedirs(num+"mask")
    except:
        pass
    os.chdir(num+"mask")

    try:
        bone = msh.Mesh(f[0])
        bone.verts[:,3] = 1
        bone.fondre(msh.Mesh(f[1]))
        bone.write("toto.mesh")
        err = os.system(exe.tetgen + " -pgaq1.2ANEF toto.mesh  > /dev/null 2>&1")
        if err:
            print "Error with tetgen on case " + num
            try:
                bone.writeOBJ("toto.obj")

        mesh = msh.Mesh("toto.1.mesh")
        ext_point_ind = mesh.tris[mesh.tris[:,-1]==0][0][0]
        ext_ref=None
        for t in mesh.tets:
            if ext_point_ind in t:
                ext_ref = t[-1]
        mesh.tris = mesh.tris[mesh.tris[:,-1]>0]
        mesh.tets = mesh.tets[mesh.tets[:,-1]==ext_ref]
        mesh.discardUnused()
        mesh.write(os.path.join(os.path.dirname(f[0]), num+".mask.mesh"))
        os.chdir("..")
    except:
        print "Error while generating mask for case " + num
        os.chdir("..")
    shutil.rmtree(num+"mask")


if __name__=="__main__":

    #Folders and static files
    csvFile        = "/home/norgeot/dev/own/FaciLe/pipeline/liste.csv"
    boneTemplate   = "/home/norgeot/dev/own/FaciLe/OsTemplate.mesh"
    massTemplate   = "/home/norgeot/dev/own/FaciLe/MassTemplate.mesh"
    mandTemplate   = "/home/norgeot/dev/own/FaciLe/MandTemplate.mesh"
    templateSphere = "/home/norgeot/dev/own/FaciLe/projects/warping/demo/sphere.o1.mesh"
    mshScript      = "/home/norgeot/dev/own/FaciLe/pipeline/cleanSTL.mlx"

    #Output folders
    outFolder      = "/Data/Facile2"
    remeshedFolder = os.path.join(outFolder,"remeshed")
    musclesFolder  = os.path.join(outFolder, "muscles")
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    if not os.path.exists(remeshedFolder):
        os.makedirs(remeshedFolder)
    if not os.path.exists(musclesFolder):
        os.makedirs(musclesFolder)

    #FTP info
    ftpUsr   = raw_input("Enter your ftp username:")
    ftpPwd   = getpass.getpass()
    IPadress = "134.157.66.224"
    ftpDir   = "Projets/FaciLe/Data/AllDataRaw"

    #Get the files present on the FTP directory if they do not exist yet
    ftp = FTP(IPadress, ftpUsr, ftpPwd)
    ftp.cwd(ftpDir)
    files = [ f for f in ftp.nlst() if ".mesh" in f or ".stl" in f ]
    files = [ f for f in files if not os.path.exists(os.path.join(outFolder,f.split(".")[1].zfill(3) + "_" + newName(f) + ".mesh")) ]
    files.sort(key=lambda f: int(f.split(".")[1]))
    print str(len(files)) + " files to copy from " + IPadress + ":" + ftpDir
    if len(files)>0:
        run(copy,files)

    #Convert the files to a .mesh format
    files = [os.path.join(outFolder, f) for f in os.listdir(outFolder) if ".stl" in f]
    files.sort()
    print str(len(files)) + " files to convert from .stl to .mesh"
    if len(files)>0:
        run(convertToMesh, files)

    #Remove duplicates and discard unused from mesh files
    if len(files)>0:
        files = [os.path.join(outFolder, f) for f in os.listdir(outFolder) if ".mesh" in f]
        files.sort()
        print str(len(files)) + " files to clean"
        if len(files)>0:
            run(cleanMesh, files)

    #Files to rescale and merge together for the bone structure
    files = [f for f in os.listdir(outFolder) if ".mesh" in f]
    files = [f for f in files if "015" not in f]
    files.sort()
    groups = group(files, condition2)
    groups = [g for g in groups if not os.path.exists(os.path.join(remeshedFolder, g[0][:3] +"_bone.mesh")) ]
    print str(len(groups)) + " groups to scale and merge"
    if len(groups)>0:
        run(process, groups)

    #Files to remesh
    files = [ os.path.join(remeshedFolder, f) for f in os.listdir(remeshedFolder) if len(f.split("."))==2 and ".mesh" in f ]
    files = [ f for f in files if not os.path.exists(f[:-5]+".o.mesh") ]
    good = readCSV(csvFile)
    good[67]   = 0
    files = [f for f in files if good[int(f.split("/")[-1][:3])-1]]
    files.sort()
    print str(len(files)) + " files to remesh"
    if len(files)>0:
        run(remesh, files)

    #Files to align
    files = [ f for f in os.listdir(remeshedFolder) if ".o.mesh" in f and not os.path.exists(os.path.join(remeshedFolder,f.split(".")[0] + ".o.final.mesh")) ]
    files.sort()
    groups = [ [ os.path.join(remeshedFolder, f) for f in files if int(f[:3])==i ] for i in range(100) ]
    groups = [g for g in groups if len(g)]
    print str(len(groups)) + " groups to align"
    if len(groups)>0:
        run(align, groups)

    #Files to warp
    files = [ os.path.join(remeshedFolder,f) for f in os.listdir(remeshedFolder) if "bone.o.final.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1].split(".")[0] + ".warped.mesh" not in os.listdir(remeshedFolder) ]
    files.sort()
    print str(len(files)) + " files to warp"
    if len(files)>0:
        run(warp, files)

    #Files to compute the signedDistance
    files = [ os.path.join(remeshedFolder,f) for f in os.listdir(remeshedFolder) if "bone.warped.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1][:3] + "_bone.signed.mesh" not in os.listdir(remeshedFolder) ]
    files.sort()
    print str(len(files)) + " files to compute the signed distance on"
    if len(files)>0:
        run(signedDistance, files)

    #Generate the masks
    faces = [ os.path.join(remeshedFolder,f) for f in os.listdir(remeshedFolder) if "face.o.final.mesh" in f]
    bones = [ os.path.join(remeshedFolder,f) for f in os.listdir(remeshedFolder) if "bone.warped.mesh"  in f]
    faces.sort()
    bones.sort()
    files = [[b,f] for b,f in zip(bones, faces)]
    files.sort()
    print str(len(files)) + " files to generate a mask from"
    if len(files)>0:
        run(generateMask, files)


    """
    #MASSETERS AND MANDIBULES
    massFiles = [f for f in os.listdir(outFolder) if "mass.mesh" in f]
    mandFiles = [f for f in os.listdir(outFolder) if "mand.mesh" in f]
    massFiles.sort()
    mandFiles.sort()
    files = massFiles + mandFiles
    files = [f for f in files if not f[:-5] + ".R.mesh" in "".join(os.listdir(musclesFolder))]
    #Processing the masseters and mandibules
    print str(len(files)) + " files to cut in half"
    if len(files)>0:
        run(cut, files)

    align=False
    liste = os.listdir(musclesFolder)
    files = [os.path.join(musclesFolder,f) for f in liste if (".R." in f or ".L." in f)]
    files = [f for f in files  if f.split("/")[-1][:-5] + ".final." not in "".join(liste) and f.split("/")[-1][:-5] + ".warped." not in "".join(liste)]
    files.sort()

    if len(files)>0:
        if align:
            print str(len(files)) + " halved files to align"
            run(alignToTemplate,files)
            print str(len(files)) + " halved files to warp"
            run(warpWithBlender, files)

    files = [os.path.join(musclesFolder,f) for f in os.listdir(musclesFolder) if ".warped" in f]
    print str(len(files)) + " warped objects to compute the signed distance on"
    if len(files)>0:
        run(signedDistance, files)


    #liste = os.listdir(musclesFolder)
    #files = [os.path.join(musclesFolder,f) for f in liste if (".R.final.mesh" in f or ".L.final.mesh" in f) and f[:-5] + ".final.o.mesh" not in "".join(liste)]
    #print str(len(files)) + " halved files to remesh"
    #if len(files)>0:
    #    run(remesh,files)

    """
