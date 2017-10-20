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
import tempfile
from functools import wraps
import string
import random
from scipy.spatial.distance import cdist

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
def nearest_neighbor(src, dst):
    all_dists = cdist(src, dst, 'euclidean')
    indices = all_dists.argmin(axis=1)
    distances = all_dists[np.arange(all_dists.shape[0]), indices]
    return distances, indices
def adapt_box_to(f, maxNb=20000):
    shutil.copyfile("/home/norgeot/box.mesh","box.mesh")
    cube = msh.Mesh("box.mesh")
    mesh = msh.Mesh(f)
    step = 1 if len(mesh.verts)<maxNb else int(len(mesh.verts)/maxNumPoints)+1
    dists, _ = nearest_neighbor(cube.verts[:,:3], mesh.verts[::step,:3])
    cube.scalars = np.array(dists)
    cube.scaleSol(0.002, 1, absolute=True)
    cube.write("box.1.mesh")
    cube.writeSol("box.1.sol")
    err = os.system("mmg3d_O3 box.1.mesh -hgrad 1.5  > /dev/null 2>&1")
    if err:
        raise FacileError("mmg3d failure")

#Decorator
class FacileError(Exception):
    pass
def debug():
    def true_decorator(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            r=None
            print f.__name__ + " : "+'\033[94m'+"RUNNING"+'\033[0m'+" on " + str(*args)
            try:
                tmpdir = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(8)])
                try:
                    os.makedirs(tmpdir)
                except:
                    print "ERROR generating temporary directory"
                    return None
                os.chdir(tmpdir)
                r = f(*args, **kwargs)
            except Exception as e:
                print f.__name__ + " : "+'\033[91m'+"FAILURE"+'\033[0m'+" on " + str(*args) + ": " + type(e).__name__ + ": " + str(e)
                pass
            else:
                print f.__name__ + " : "+'\033[92m'+"SUCCESS"+'\033[0m'+" on " + str(*args)
            finally:
                os.chdir("..")
                shutil.rmtree(tmpdir)
                return r
        return wrapped
    return true_decorator
dryRun  = False # Run just for fun
oneStep = False # If only wish to process one step
def run(func, liste, parallel=True, maxi=64):
    if len(liste)>0:
        num = min( maxi, min(len(liste), mp.cpu_count()-1 )) if parallel else 1
        print '\033[95m' + "## EXECUTING '" + func.__name__ + "' on " + str(len(liste)) + " cases on " + str(num) + " process(es)." + '\033[0m'
        if not dryRun:
            res = []
            if parallel:
                pool = mp.Pool(processes=num)
                res = pool.map(func, liste )
            else:
                res = []
                for l in liste:
                    res.append(func(l))
            if oneStep:
                print '\033[95m' + "ONESTEP -> EXITING..." + '\033[0m'
                sys.exit()
            return res
        else:
            return 0
    else:
        print '\033[95m' + "## 0 arguments to execute '" + func.__name__ + "' on, skipping..." + '\033[0m'

#Functions
@debug()
def ftpCopy(f):
    ftp = FTP(IPadress, ftpUsr, ftpPwd)
    ftp.cwd(ftpDir)
    num = f.split(".")[1].zfill(3)
    localFile = os.path.join(rawFolder, num + "_" + newName(f) + "." + f.split(".")[-1])
    if not num + "_" + newName(f) + ".mesh" in "".join(os.listdir(rawFolder)):
        if not os.path.isfile(localFile):
            with open(localFile, 'wb') as ff:
                ftp.retrbinary('RETR %s' % f, ff.write)
@debug()
def convertToMesh(f):
    if "stl" in f:
        if not stlToMesh(f):
            raise FacileError("conversion failure")
@debug()
def cleanMesh(f):
    mesh = msh.Mesh()
    mesh.get_infos(f)
    nV = mesh.numItems[0]
    nT = mesh.numItems[1]
    if nV>nT:
        mesh = msh.Mesh(f)
        mesh.discardDuplicateVertices()
        mesh.discardUnused()
        mesh.write(f)
@debug()
def process(g):
    newBone = os.path.join(remeshedFolder, g[0][:3] + "_bone.mesh")
    newFace = os.path.join(remeshedFolder, g[0][:3] + "_face.mesh")
    newMass = os.path.join(remeshedFolder, g[0][:3] + "_mass.mesh")

    faceFile = [f for f in g if "face" in f][0]
    face = msh.Mesh(os.path.join(rawFolder, faceFile))
    boneFiles = [f for f in g if f!=faceFile and "mass" not in f]
    bone   = None
    center = [0,0,0]
    bone = msh.Mesh(os.path.join(rawFolder,boneFiles[0]))
    if len(boneFiles)>1:
        for f in boneFiles[1:]:
            bone.fondre(msh.Mesh(os.path.join(rawFolder,f)))
    center   = bone.center
    scale    = 0.0035
    meshes = [face, bone]
    mass = None
    if "mass" in "".join(g):
        mass = msh.Mesh(os.path.join(rawFolder, [f for f in g if "mass" in f][0]))
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
@debug()
def remesh(f, hausd=0.005):
    err = os.system(exe.mmgs + " " + f + " -nr -nreg -hausd " + str(hausd) + " -o " + f[:-5] + ".o.mesh > /dev/null 2>&1")
    if err:
        raise FacileError("mmgs failure")
    return 0
@debug()
def align(g):
    boneFile = [f for f in g if "bone" in f][0]
    faceFile = [f for f in g if "face" in f][0]
    massFile = [f for f in g if "mass" in f][0] if "mass" in "".join(g) else ""
    num = boneFile.split("/")[-1][:3]

    err = os.system(exe.align + " -i " + boneFile + " " + boneTemplate + " -d 0.1 -o 0.95  > "+num+".txt")#/dev/null 2>&1")
    if err:
        raise FacileError("alignement failure")

    bone = msh.Mesh(boneFile)
    bone.applyMatrix(matFile = "mat_Super4PCS.txt")
    bone.applyMatrix(matFile = "mat_ICP.txt")
    bone.write(os.path.join(alignedFolder, num+"_bone.aligned.mesh"))

    err = os.system(exe.pythonICP + " -s " + os.path.join(alignedFolder, num+"_bone.aligned.mesh") + " -t " + boneTemplate + " -m mat_pyICP.txt >> " + num + ".txt")
    if err:
        pass # Cannot run python alignement...
    else:
        bone.applyMatrix(matFile = "mat_pyICP.txt")
        bone.write(os.path.join(alignedFolder, num+"_bone.aligned.mesh"))

    face = msh.Mesh(faceFile)
    face.applyMatrix(matFile = "mat_Super4PCS.txt")
    face.applyMatrix(matFile = "mat_ICP.txt")
    if not err:
        face.applyMatrix(matFile = "mat_pyICP.txt")
    face.write(os.path.join(alignedFolder, num+"_face.aligned.mesh"))
    if ".mesh" in massFile:
        mass = msh.Mesh(massFile)
        mass.applyMatrix(matFile = "mat_Super4PCS.txt")
        mass.applyMatrix(matFile = "mat_ICP.txt")
        if not err:
            mass.applyMatrix(matFile = "mat_pyICP.txt")
        mass.write(os.path.join(alignedFolder, num+"_mass.aligned.mesh"))
    return 0
@debug()
def warp(f):
    num = f.split("/")[-1][:3]

    os.system("cp " + templateSphere + " ./sphere.mesh")
    err = os.system( exe.warping + " " + f + " -p -load 40 > warping.txt" )
    if err:
        try:
            mesh = msh.Mesh(f)
            mesh.tris = np.array([ t for t in mesh.tris if np.max(t)<=len(mesh.verts) ])
            mesh.discardUnused()
            mesh.write(f)
            err = os.system( exe.warping + " " + f + " -p -load 70 > warping.txt" )
            if err:
                raise FacileError("cleaned warping failure")
        except:
            raise FacileError("warping failure")

    warped = msh.Mesh("sphere.d.mesh")
    ext_ref = 2
    warped.tris = warped.tris[warped.tris[:,-1] != ext_ref]
    warped.tets = np.array([])
    warped.discardUnused()
    warped.write(os.path.join(warpedFolder,f.split("/")[-1].split(".")[0] + ".warped.mesh"))
    return 0
@debug()
def signedDistance(f):
    adapt = True
    if adapt:
        adapt_box_to(f)
    else:
        cube=msh.Mesh(cube=[0,1,0,1,0,1])
        cube.write("box.mesh")
        err = os.system( exe.tetgen + " -pgANEF box.mesh > /dev/null 2>&1")
        if err:
            raise FacileError('tetgen failure')
        err = os.system( exe.mmg3d + " box.1.mesh -hausd 0.04 -hmax 0.04 > /dev/null 2>&1" )
        if err:
            raise FacileError('mmg3d failure')

    err = os.system( exe.mshdist + " -ncpu 16 -noscale box.1.o.mesh " + f + " > /dev/null 2>&1")
    if err:
        raise FacileError('mshdist failure')
    name = f.split("/")[-1].split(".")[0]
    os.system("mv box.1.o.mesh " + os.path.join(signedFolder, name + ".signed.mesh"))
    os.system("mv box.1.o.sol  " + os.path.join(signedFolder, name + ".signed.sol"))
    return 0
@debug()
def cut(f):
    mesh=msh.Mesh(os.path.join(rawFolder, f))
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

def alignToTemplate(f):
    num = f.split("/")[-1][:3]
    err=0
    """
    if "mass" in f:
        err = os.system(exe.align + " -i " + f + " " + massTemplate + " -d 0.1 -o 0.95  > "+num+"_mass.txt")#/dev/null 2>&1")
    elif "mand" in f:
        err = os.system(exe.align + " -i " + f + " " + mandTemplate + " -d 0.1 -o 0.95  > "+num+"_mand.txt")#/dev/null 2>&1")
    if err:
        print "-- Error aligning " + boneFile
    """
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
    return 0

@debug()
def warpWithBlender(f,resolution=41):
    err1 = os.system( exe.boundingMesh + " -i " + f + " -o carved.mesh -r " + str(resolution) + " > /dev/null 2>&1")
    err2 = os.system("blender --background --python blender_warp.py -- carved.mesh " + f + " warped.mesh > /dev/null 2>&1")
    err3 = os.system(exe.mmgs + " warped.mesh -o warped.o.mesh -hausd 0.002 > /dev/null 2>&1")
    if (err1+err2+err3):
        raise FacileError("blender warping failure")
    shutil.copyfile("warped.o.mesh", f[:-5]+".rewarped.mesh")

def cleanWithMeshlab(f, new=None):
    mesh = msh.Mesh(f)
    mesh.writeOBJ("tmp.obj")
    err = os.system("LC_ALL=C meshlabserver -i tmp.obj -o cleaned.obj -s " + intScript + " > /dev/null 2>&1" )
    if err:
        raise FacileError("meshlab failure")
    mesh = obj2Mesh("cleaned.obj")
    if new is not None:
        mesh.write(new)
    else:
        mesh.write(f)
@debug()
def intersects(f):
    os.system(exe.tetgen + " -d " + f + " > log.txt")
    with open("log.txt","r") as f:
        if "No faces are intersecting" in "".join(f.readlines()):
            return False
        else:
            return True
@debug()
def generateMask(f):
    num = f[0].split("/")[-1][:3]

    bone, face = None, None

    if intersects(f[0]):
        cleanWithMeshlab(f[0],new="bone.mesh")
        if intersects("bone.mesh"):
            raise FacileError("Bone intersect")
        bone = msh.Mesh("bone.mesh")
    else:
        bone = msh.Mesh(f[0])

    if intersects(f[1]):
        cleanWithMeshlab(f[1],new="face.mesh")
        if intersects("face.mesh"):
            raise FacileError("Face intersect")
        face = msh.Mesh("face.mesh")
    else:
        face = msh.Mesh(f[1])

    bone.tris[:,3] = 1
    face.tris[:,3] = 2
    face.discardUnused()
    bone.fondre(face)
    bone.write("mask.mesh")

    if intersects("mask.mesh"):
        cleanWithMeshlab("mask.mesh")
        if intersects("mask.mesh"):
            raise FacileError('Bone and Face intersect')

    err = os.system(exe.tetgen + " -pgaAYNEF mask.mesh > /dev/null 2>&1")
    if err:
        raise FacileError("tetgen error")

    mesh = msh.Mesh("mask.1.mesh")

    #Exterior point = closest to [0,0,0]
    ext_point_ind = np.argmin([np.linalg.norm(x) for x in mesh.verts[:,:3]])
    ext_ref=None
    for t in mesh.tets:
        if ext_point_ind in t[:4]:
            ext_ref = t[-1]
            break
    mesh.tets = mesh.tets[mesh.tets[:,-1]==ext_ref]
    mesh.tets[:,4] = 2

    for t in mesh.tris:
        if ext_point_ind in t[:3]:
            ext_ref = t[-1]
            break
    mesh.tris = mesh.tris[mesh.tris[:,3]>0]
    M = mesh.tris[:,-1]==ext_ref
    mesh.tris[M==1][:,3] = 1
    mesh.tris[M==0][:,3] = 0

    mesh.discardUnused()

    mesh.write(os.path.join(maskedFolder, num + "_mask.mesh"))

if __name__=="__main__":

    #Folders and static files
    csvFile        = "/home/norgeot/dev/own/FaciLe/pipeline/liste.csv"
    boneTemplate   = "/home/norgeot/dev/own/FaciLe/OsTemplate.mesh"
    massTemplate   = "/home/norgeot/dev/own/FaciLe/MassTemplate.mesh"
    mandTemplate   = "/home/norgeot/dev/own/FaciLe/MandTemplate.mesh"
    templateSphere = "/home/norgeot/dev/own/FaciLe/projects/warping/demo/sphere.o1.mesh"
    mshScript      = "/home/norgeot/dev/own/FaciLe/pipeline/cleanSTL.mlx"
    intScript      = "/home/norgeot/dev/own/FaciLe/pipeline/cleanIntersections.mlx"

    #Output folders
    outFolder      = "/Data/Facile3"
    rawFolder      = os.path.join(outFolder, "raw")
    remeshedFolder = os.path.join(outFolder, "remeshed")
    alignedFolder  = os.path.join(outFolder, "aligned")
    warpedFolder   = os.path.join(outFolder, "warped")
    signedFolder   = os.path.join(outFolder, "signed")
    maskedFolder   = os.path.join(outFolder, "masked")
    morphedFolder  = os.path.join(outFolder, "morphed")
    musclesFolder  = os.path.join(outFolder, "muscles")
    if not os.path.exists(rawFolder):
        os.makedirs(rawFolder)
    if not os.path.exists(remeshedFolder):
        os.makedirs(remeshedFolder)
    if not os.path.exists(musclesFolder):
        os.makedirs(musclesFolder)
    if not os.path.exists(alignedFolder):
        os.makedirs(alignedFolder)
    if not os.path.exists(warpedFolder):
        os.makedirs(warpedFolder)
    if not os.path.exists(signedFolder):
        os.makedirs(signedFolder)
    if not os.path.exists(maskedFolder):
        os.makedirs(maskedFolder)
    if not os.path.exists(morphedFolder):
        os.makedirs(morphedFolder)

    #FTP info
    ftpUsr   = "lnorgeot"#raw_input("Enter your ftp username:")
    ftpPwd   = "Lolopolo29**"#getpass.getpass()
    IPadress = "134.157.66.224"
    ftpDir   = "Projets/FaciLe/Data/AllDataRaw"

    # 1 - Copy from ftp to rawFolder
    ftp = FTP(IPadress, ftpUsr, ftpPwd)
    ftp.cwd(ftpDir)
    f_ftp = [ f for f in ftp.nlst() if ".mesh" in f or ".stl" in f ]
    f_ftp = [ f for f in f_ftp if not os.path.exists(os.path.join(rawFolder,f.split(".")[1].zfill(3) + "_" + newName(f) + ".mesh")) ]
    f_ftp.sort(key=lambda f: int(f.split(".")[1]))
    run(ftpCopy, f_ftp)

    # 2 - Convert from .stl to .mesh in rawFolder
    files = [os.path.join(rawFolder, f) for f in os.listdir(rawFolder) if ".stl" in f]
    files.sort()
    run(convertToMesh, files)

    # 3 - Clean the meshes in rawFolder
    files = [os.path.join(rawFolder, f) for f in os.listdir(rawFolder) if ".mesh" in f] if len(files) else []
    files.sort()
    run(cleanMesh, files)

    # 4 - Scale files and merge the bones to remeshedFolder
    files = [f for f in os.listdir(rawFolder) if ".mesh" in f]
    files = [f for f in files if "015" not in f]
    files.sort()
    groups = group(files, condition2)
    groups = [g for g in groups if not os.path.exists(os.path.join(remeshedFolder, g[0][:3] +"_bone.mesh")) ]
    run(process, groups)

    # 5 - Remesh the files
    files = [ os.path.join(remeshedFolder, f) for f in os.listdir(remeshedFolder) if len(f.split("."))==2 and ".mesh" in f ]
    files = [ f for f in files if not os.path.exists(f[:-5]+".o.mesh") ]
    good = readCSV(csvFile)
    good[67]   = 0
    files = [f for f in files if good[int(f.split("/")[-1][:3])-1]]
    files.sort()
    run(remesh, files)

    # 6 - Align files
    files = [ f for f in os.listdir(remeshedFolder) if ".o.mesh" in f and not os.path.exists(os.path.join(alignedFolder,f.split(".")[0] + ".aligned.mesh")) ]
    files.sort()
    groups = [ [ os.path.join(remeshedFolder, f) for f in files if int(f[:3])==i ] for i in range(500) ]
    groups = [g for g in groups if len(g)]
    run(align, groups)

    # 7 - Warp the bones
    files = [ os.path.join(alignedFolder,f) for f in os.listdir(alignedFolder) if "bone.aligned.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1].split(".")[0] + ".warped.mesh" not in os.listdir(warpedFolder) ]
    files.sort()
    run(warp, files)

    # 8 - Compute the signed distance on bones
    files = [ os.path.join(warpedFolder,f) for f in os.listdir(warpedFolder) if "bone.warped.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1][:3] + "_bone.signed.mesh" not in os.listdir(signedFolder) ]
    files.sort()
    run(signedDistance, files, maxi=4)

    # 9 - Compute the signed distance on faces
    """
    files = [ os.path.join(alignedFolder,f) for f in os.listdir(alignedFolder) if "face.aligned.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1][:3] + "_face.signed.mesh" not in os.listdir(signedFolder) ]
    files.sort()
    run(signedDistance, files, maxi=4)
    """

    # 10 - Generate the masks
    faces = [ os.path.join(alignedFolder,f) for f in os.listdir(alignedFolder) if "face.aligned.mesh" in f]
    bones = [ os.path.join(warpedFolder,f) for f in os.listdir(warpedFolder) if "bone.warped.mesh"  in f]
    faces.sort()
    bones.sort()
    files = [[b,f] for b,f in zip(bones, faces)]
    files.sort()
    run(generateMask, files)

    """
    f = "/Data/Facile3/raw/045_face.mesh"
    err1 = os.system( exe.boundingMesh + " -i " + f + " -o carved.mesh -r 51")
    err2 = os.system("blender --background --python blender_warp.py -- carved.mesh " + f + " warped.mesh")
    err3 = os.system(exe.mmgs + " -nr -nreg warped.mesh -o warped.o.mesh -hausd 0.002")
    newFaces = [ os.path.join(remeshedFolder,f) for f in os.listdir(remeshedFolder) if "face.o.final.rewarped.mesh" in f]
    for f in newFaces:
        print f
    """



    """
    #MASSETERS AND MANDIBULES
    massFiles = [f for f in os.listdir(rawFolder) if "mass.mesh" in f]
    mandFiles = [f for f in os.listdir(rawFolder) if "mand.mesh" in f]
    massFiles.sort()
    mandFiles.sort()
    files = massFiles + mandFiles
    files = [f for f in files if not f[:-5] + ".R.mesh" in "".join(os.listdir(musclesFolder))]
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
    run(signedDistance, files)

    #liste = os.listdir(musclesFolder)
    #files = [os.path.join(musclesFolder,f) for f in liste if (".R.final.mesh" in f or ".L.final.mesh" in f) and f[:-5] + ".final.o.mesh" not in "".join(liste)]
    #print str(len(files)) + " halved files to remesh"
    #if len(files)>0:
    #    run(remesh,files)
    """
