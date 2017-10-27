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
import time

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
    shutil.copyfile(templates["box"],"box.mesh")
    cube = msh.Mesh("box.mesh")
    mesh = msh.Mesh(f)
    step = 1 if len(mesh.verts)<maxNb else int(len(mesh.verts)/maxNb)+1
    dists, _ = nearest_neighbor(cube.verts[:,:3], mesh.verts[::step,:3])
    cube.scalars = np.array(dists)
    cube.scaleSol(0.001, 0.5, absolute=True)
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
            t = time.time()
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
                print f.__name__ + " : "+'\033[91m'+"FAILURE"+'\033[0m'+" on " + str(*args) + ": " + type(e).__name__ + ": " + str(e) + ", in " + str(int(time.time() - t)) + " s"
                pass
            else:
                print f.__name__ + " : "+'\033[92m'+"SUCCESS"+'\033[0m'+" on " + str(*args) + ", in " + str(int(time.time() - t)) + " s"
            finally:
                os.chdir("..")
                shutil.rmtree(tmpdir)
                return r
        return wrapped
    return true_decorator
dryRun  = True  # Run just for fun
oneStep = False  # If only wish to process one step
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
    localFile = os.path.join(directories["raw"], num + "_" + newName(f) + "." + f.split(".")[-1])
    if not num + "_" + newName(f) + ".mesh" in "".join(os.listdir(directories["raw"])):
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
def scale(g):
    faceFile  = [f for f in g if "face" in f][0]
    boneFiles = [f for f in g if f!=faceFile and "mass" not in f]
    bone = msh.Mesh(os.path.join(directories["raw"],boneFiles[0]))
    if len(boneFiles)>1:
        for f in boneFiles[1:]:
            bone.fondre(msh.Mesh(os.path.join(directories["raw"],f)))
    center   = bone.center
    scale    = 0.0035
    for f in g:
        mesh = msh.Mesh(os.path.join(directories["raw"],f))
        mesh.verts[:,:3] -= center
        mesh.verts[:,:3] *= scale
        mesh.verts[:,:3] += [0.5,0.5,0.5]
        mesh.write(os.path.join(directories["scaled"],f))
@debug()
def remesh(f, hausd=0.0025):
    err = os.system(exe.mmgs + " " + os.path.join(directories["scaled"],f) + " -nr -nreg -hausd " + str(hausd) + " -o " + os.path.join(directories["remeshed"],f) + " > /dev/null 2>&1")
    if err:
        raise FacileError("mmgs failure")
    return 0
@debug()
def merge(g):
    newBone = os.path.join(directories["merged"], g[0][:3] + "_bone.mesh")
    newFace = os.path.join(directories["merged"], g[0][:3] + "_face.mesh")
    newMass = os.path.join(directories["merged"], g[0][:3] + "_mass.mesh")

    faceFile = [f for f in g if "face" in f][0]
    face = msh.Mesh(os.path.join(directories["remeshed"], faceFile))
    boneFiles = [f for f in g if f!=faceFile and "mass" not in f]

    if len(boneFiles)==0:
        return 1

    bone = msh.Mesh(os.path.join(directories["remeshed"],boneFiles[0]))
    if len(boneFiles)>1:
        for f in boneFiles[1:]:
            bone.fondre(msh.Mesh(os.path.join(directories["remeshed"],f)))

    if "mass" in "".join(g):
        mass = msh.Mesh(os.path.join(directories["remeshed"], [f for f in g if "mass" in f][0]))
        mass.write(newMass)
    bone.write(newBone)
    face.write(newFace)
@debug()
def align(g):
    boneFile = [f for f in g if "bone" in f][0]
    faceFile = [f for f in g if "face" in f][0]
    massFile = [f for f in g if "mass" in f][0] if "mass" in "".join(g) else ""
    num = boneFile.split("/")[-1][:3]

    err = os.system(exe.align + " -i " + boneFile + " " + templates["bone"] + " -d 0.1 -o 0.95  > "+num+".txt")#/dev/null 2>&1")
    if err:
        raise FacileError("alignement failure")

    bone = msh.Mesh(boneFile)
    bone.applyMatrix(matFile = "mat_Super4PCS.txt")
    bone.applyMatrix(matFile = "mat_ICP.txt")
    bone.write(os.path.join(directories["aligned"], num+"_bone.mesh"))

    err = os.system(exe.pythonICP + " -s " + os.path.join(directories["aligned"], num+"_bone.mesh") + " -t " + templates["bone"] + " -m mat_pyICP.txt >> " + num + ".txt")
    if err:
        pass # Cannot run python alignement...
    else:
        bone.applyMatrix(matFile = "mat_pyICP.txt")
        bone.write(os.path.join(directories["aligned"], num+"_bone.mesh"))

    face = msh.Mesh(faceFile)
    face.applyMatrix(matFile = "mat_Super4PCS.txt")
    face.applyMatrix(matFile = "mat_ICP.txt")
    if not err:
        face.applyMatrix(matFile = "mat_pyICP.txt")
    face.write(os.path.join(directories["aligned"], num+"_face.mesh"))
    if ".mesh" in massFile:
        mass = msh.Mesh(massFile)
        mass.applyMatrix(matFile = "mat_Super4PCS.txt")
        mass.applyMatrix(matFile = "mat_ICP.txt")
        if not err:
            mass.applyMatrix(matFile = "mat_pyICP.txt")
        mass.write(os.path.join(directories["aligned"], num+"_mass.mesh"))
    return 0
@debug()
def warp(f):
    num = f.split("/")[-1][:3]

    os.system("cp " + templates["sphere"] + " ./sphere.mesh")
    err = os.system( exe.warping + " " + f + " -p -nit 150 -load 40 > warping.txt" )
    if err:
        raise FacileError("Warping failure")
    """
        try:
            mesh = msh.Mesh(f)
            mesh.tris = np.array([ t for t in mesh.tris if np.max(t)<=len(mesh.verts) ])
            mesh.discardUnused()
            mesh.write(f)
            err = os.system( exe.warping + " " + f + " -p -nit 150 -load 40 > warping.txt" )
            if err:
                raise FacileError("cleaned warping failure")
        except:
            raise FacileError("warping failure")
    """
    warped = msh.Mesh("sphere.d.mesh")
    ext_ref = 2
    warped.tris = warped.tris[warped.tris[:,-1] != ext_ref]
    warped.tets = np.array([])
    warped.discardUnused()
    warped.write(os.path.join(directories["warped"],f.split("/")[-1].split(".")[0] + ".warped.mesh"))
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
    os.system("mv box.1.o.mesh " + os.path.join(directories["signed"], name + ".mesh"))
    os.system("mv box.1.o.sol  " + os.path.join(directories["signed"], name + ".sol"))
    return 0
@debug()
def cut(f):
    mesh=msh.Mesh(os.path.join(directories["raw"], f))
    MAT = mesh.toUnitMatrix()
    mesh.applyMatrix(MAT)
    newF = os.path.join(directories["muscles"], f[:-5] + ".scaled.mesh")
    mesh.write(newF)
    err = remesh(newF, hausd=0.001)
    mesh = msh.Mesh(os.path.join(directories["muscles"], f[:-5] + ".scaled.o.mesh"))
    half1, half2 = cutMeshInHalf(mesh)
    half1.write(os.path.join(directories["muscles"], f[:-5]+".R.mesh"))
    half2.write(os.path.join(directories["muscles"], f[:-5]+".L.mesh"))
    os.remove(os.path.join(directories["muscles"], f[:-5] + ".scaled.mesh"))
    os.remove(os.path.join(directories["muscles"], f[:-5] + ".scaled.o.mesh"))
    os.remove(os.path.join(directories["muscles"], f[:-5] + ".scaled.o.sol"))
def alignToTemplate(f):
    num = f.split("/")[-1][:3]
    err=0
    """
    if "mass" in f:
        err = os.system(exe.align + " -i " + f + " " + templates["masseter"] + " -d 0.1 -o 0.95  > "+num+"_mass.txt")#/dev/null 2>&1")
    elif "mand" in f:
        err = os.system(exe.align + " -i " + f + " " + templates["mandible"] + " -d 0.1 -o 0.95  > "+num+"_mand.txt")#/dev/null 2>&1")
    if err:
        print "-- Error aligning " + boneFile
    """
    mesh = msh.Mesh(f)
    #mesh.applyMatrix(matFile = "mat_Super4PCS.txt")
    #mesh.applyMatrix(matFile = "mat_ICP.txt")
    #mesh.write(f[:-5]+".aligned.mesh")
    todo = f#f[:-5]+".aligned.mesh"
    if "mass" in f:
        err = os.system(exe.pythonICP + " -s " + todo + " -t " + templates["masseter"] + " -m mat_pyICP.txt >> " + num + "_mass.txt")
    elif "mand" in f:
        err = os.system(exe.pythonICP + " -s " + todo + " -t " + templates["mandible"] + " -m mat_pyICP.txt >> " + num + "_mand.txt")
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
def intersects(f):
    err = os.system(exe.tetgen + " -d " + f + " > log.txt")
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
        #cleanWithMeshlab(f[0],new="bone.mesh")
        #if intersects("bone.mesh"):
        raise FacileError("Bone intersect")
        #bone = msh.Mesh("bone.mesh")
    else:
        bone = msh.Mesh(f[0])

    if intersects(f[1]):
        #cleanWithMeshlab(f[1],new="face.mesh")
        #if intersects("face.mesh"):
        raise FacileError("Face intersect")
        #face = msh.Mesh("face.mesh")
    else:
        face = msh.Mesh(f[1])

    face.tris = face.tris[face.tris[:,-1]!=2]
    face.tris[:,3] = 2
    face.discardUnused()
    bone.tris[:,3]  = 1
    bone.fondre(face)
    bone.write("mask.mesh")

    if intersects("mask.mesh"):
        #cleanWithMeshlab("mask.mesh")
        #if intersects("mask.mesh"):
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

    for t in mesh.tris:
        if t[-1]==1:
            for i in t[:3]:
                mesh.verts[i,-1]=1

    mesh.write(os.path.join(directories["masked"], num + "_mask.mesh"))

    err = os.system(exe.mmg3d + " " + os.path.join(directories["masked"], num + "_mask.mesh") + " -o " + os.path.join(directories["masked"], num + "_mask.o.mesh") + " -hausd 0.0005 -nosurf -hgrad 1.15 > /dev/null 2>&1")
    if err:
        raise FacileError("mmg3d error")
    os.system("rm " + os.path.join(directories["masked"], num + "_mask.o.sol"))


    template = msh.Mesh(templates["morphingSkull"])
    template.tets = np.array([])
    template.tris = template.tris[template.tris[:,-1]==1]
    template.discardUnused()

    mesh = msh.Mesh(os.path.join(directories["masked"], num + "_mask.o.mesh"))
    n = len(mesh.verts)
    mesh.tris = mesh.tris[mesh.tris[:,-1]==1]
    mesh.tets = np.array([])
    mesh.discardUnused()
    mesh.vectors = np.zeros((n,3))
    mesh.vectors[:len(template.verts)] = mesh.verts[:,:3] - template.verts[:,:3] # A remplacer par le resultat du morphing
    mesh.writeSol( os.path.join(directories["masked"], num + "_mask.o.sol") )
@debug()
def morph(g, nit=400):
    """
    g[0] = path to the signed distance
    g[1] = path to the template mesh
    """
    signedDist = g[0]
    templateMesh = g[1]
    os.system("cp " + templateMesh + " template.mesh")
    dRef  = [2] #Fixed surface inside the template
    elRef = [2] #Elements inside the fixed surface
    bRef  = [] #Follower elements
    cmd = " ".join((
        "morphing",
        " -dref "  + str(len(dRef))  + " " + " ".join([str(d) for d in dRef]),
        " -elref " + str(len(elRef)) + " " + " ".join([str(d) for d in elRef]),
        #" -bref "  + str(len(bRef))  + " " + " ".join([str(d) for d in bRef]),
        " -nit " + str(nit),
        " " + signedDist,
        " template.mesh",
        " > /dev/null 2>&1"
    ))
    #print cmd
    if True:
        err = os.system(cmd)
        if err:
            raise FacileError("morphing failure")
        name = signedDist.split("/")[-1].split(".")[0]
        newMesh = os.path.join(directories["morphed"], name+".mesh")
        newSol  = os.path.join(directories["morphed"], name+".sol")
        os.system("mv " + signedDist[:-5] + ".1.mesh "     + newMesh)
        os.system("mv " + signedDist[:-5] + ".1.depl.sol " + newSol)
        mesh = msh.Mesh(newMesh)
        mesh.readSol(newSol)
        mesh.tets = np.array([])
        mesh.discardUnused()
        mesh.write(newMesh)
        mesh.writeSol(newSol)








import numpy as np
import time
def scalar(d1,d2):
    return np.sum(np.multiply(d1,d2))
def read_data(filename):
    with open(filename,"r") as f:
        return np.array([[float(x) for x in l.split()] for l in f.readlines()[:4859]])
    return None
def center_data(d):
    return d - np.mean(d,axis=0)
def cov(d):
    return np.array([[scalar(x,y) for x in d] for y in d])
def eig(c):
    eVal, eVec = np.linalg.eig(c)
    eVec = eVec.transpose()
    idx = eVal.argsort()[::-1]
    eVal = eVal[idx]
    eVec = eVec[idx,:]
    return eVal, eVec
def get_principal_components(v, d):
    pc  = np.array([np.sum([v[i,j]*d[j] for j in range(len(d))],axis=0) for i in range(len(d))])
    pcn = [x/np.sqrt(scalar(x,x)) for x in pc]
    return pc, pcn
def reconstruct(pcn, d, n=None):
    alpha = np.array([[scalar(x,y) for y in pcn] for x in d])
    if n:
        return np.array([np.sum([alpha[i,j] * pcn[j] for j in range(n)], axis=0) for i in range(len(d))])
    else:
        return np.array([np.sum([alpha[i,j] * pcn[j] for j in range(len(d))], axis=0) for i in range(len(d))])
def PCA(d, u, n, debug=False):
    if debug:
        t0 = time.time()
    A = cov(d)
    if debug:
        print "Calcul de A: " + str(time.time() - t0)
        t0 = time.time()
    eVal, eVec = eig(A)
    if debug:
        print "Calcul des valeurs/vecteurs propres: " + str(time.time() - t0)
        t0 = time.time()
    PC,PCN = get_principal_components(eVec, d)
    if debug:
        print "Calcul des composantes principales: " + str(time.time() - t0)
        t0 = time.time()
    #N = reconstruct(PCN, d)
    N = reconstruct(PCN, np.append(d,[u],axis=0), n)
    if debug:
        print "Calcul de la reconstruction: " + str(time.time() - t0)
        t0 = time.time()
    return N[-1]








if __name__=="__main__":

    #Folders and static files
    csvFile        = "/home/norgeot/dev/own/FaciLe/pipeline/liste.csv"
    mshScript      = "/home/norgeot/dev/own/FaciLe/pipeline/cleanSTL.mlx"
    intScript      = "/home/norgeot/dev/own/FaciLe/pipeline/cleanIntersections.mlx"

    #templateFiles
    templates = {}
    templates["masseter"] = "/home/norgeot/dev/own/FaciLe/MassTemplate.mesh"
    templates["mandible"] = "/home/norgeot/dev/own/FaciLe/MandTemplate.mesh"
    templates["bone"]     = "/home/norgeot/dev/own/FaciLe/OsTemplate.mesh"
    templates["sphere"]   = "/home/norgeot/dev/own/FaciLe/projects/warping/demo/sphere.o1.mesh"
    templates["morphing"] = "/home/norgeot/templateMorphing.mesh"
    templates["morphingSkull"] = "/home/norgeot/templateMorphingSkull.mesh"
    templates["box"]      = "/home/norgeot/box.mesh"

    #Output folders
    outFolder      = "/Data/Facile3"
    dirNames       = [
        "raw",
        "scaled",
        "remeshed",
        "merged",
        "aligned",
        "warped",
        "signed",
        "masked",
        "morphed",
        "muscles",
        "reconstruction"
    ]
    directories    = {}
    for d in dirNames:
        directories[d] = os.path.join(outFolder, d)
        if not os.path.exists(directories[d]):
            os.makedirs(directories[d])

    #FTP info
    ftpUsr   = raw_input("Enter your ftp username:")
    ftpPwd   = getpass.getpass()
    IPadress = "134.157.66.224"
    ftpDir   = "Projets/FaciLe/Data/AllDataRaw"


    ############################################################################
    # I - COPY AND CLEAN
    ############################################################################

    # I.1 - Copy from ftp to rawFolder
    ftp = FTP(IPadress, ftpUsr, ftpPwd)
    ftp.cwd(ftpDir)
    f_ftp = [ f for f in ftp.nlst() if ".mesh" in f or ".stl" in f ]
    f_ftp = [ f for f in f_ftp if not os.path.exists(os.path.join(directories["raw"],f.split(".")[1].zfill(3) + "_" + newName(f) + ".mesh")) ]
    f_ftp.sort(key=lambda f: int(f.split(".")[1]))
    run(ftpCopy, f_ftp)

    # I.2 - Convert from .stl to .mesh in rawFolder
    files = [os.path.join(directories["raw"], f) for f in os.listdir(directories["raw"]) if ".stl" in f]
    files.sort()
    run(convertToMesh, files)

    # I.3 - Clean the meshes in rawFolder
    files = [os.path.join(directories["raw"], f) for f in os.listdir(directories["raw"]) if ".mesh" in f] if len(files) else []
    files.sort()
    run(cleanMesh, files)


    ############################################################################
    # II - UNIFORMIZE (SCALE, REMESH, MERGE AND ALIGN)
    ############################################################################

    # II.1 - Scale the files
    files = [ f for f in os.listdir(directories["raw"]) if ".mesh" in f ]
    files.sort()
    groups = group(files, condition2)
    groups = [g for g in groups if not os.path.exists(os.path.join(directories["scaled"], g[0][:3] +"_bone.mesh")) ]
    run(scale, groups)

    # II.2 - Remesh the files
    files = [ f for f in os.listdir(directories["scaled"]) if ".mesh" in f ]
    files = [ f for f in files if not os.path.exists(os.path.join(directories["remeshed"],f)) ]
    good = readCSV(csvFile)
    files = [f for f in files if good[int(f[:3])-1]]
    files.sort()
    run(remesh, files)

    # II.3 - merge the bones together
    files = [f for f in os.listdir(directories["remeshed"]) if ".mesh" in f]
    files.sort()
    groups = group(files, condition2)
    groups = [g for g in groups if not os.path.exists(os.path.join(directories["merged"], g[0][:3] +"_bone.mesh")) ]
    run(merge, groups)

    # II.4 - Align files
    files = [ f for f in os.listdir(directories["merged"]) if ".mesh" in f and not os.path.exists(os.path.join(directories["aligned"],f.split(".")[0] + ".mesh")) ]
    files.sort()
    groups = [ [ os.path.join(directories["merged"], f) for f in files if int(f[:3])==i ] for i in range(500) ]
    groups = [g for g in groups if len(g)]
    run(align, groups)


    ############################################################################
    # III - BONES
    ############################################################################

    # III.1 - Warp the bones
    files = [ os.path.join(directories["aligned"],f) for f in os.listdir(directories["aligned"]) if "bone.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1].split(".")[0] + ".warped.mesh" not in os.listdir(directories["warped"]) ]
    files.sort()
    run(warp, files)

    # III.2 - Compute the signed distance on bones
    files = [ os.path.join(directories["warped"],f) for f in os.listdir(directories["warped"]) if "bone.warped.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1][:3] + "_bone.mesh" not in os.listdir(directories["signed"]) ]
    files.sort()
    run(signedDistance, files, maxi=4)

    # III.3 - Morphing from the template skull to every skull
    files = [ os.path.join(directories["signed"],f) for f in os.listdir(directories["signed"]) if "bone.mesh" in f ]
    files = [ [f, templates["morphingSkull"]] for f in files if f.split("/")[-1][:3] + "_bone.mesh" not in os.listdir(directories["morphed"]) ]
    files.sort()
    run(morph, files)


    ############################################################################
    # IV - FACES
    ############################################################################

    # IV.1 - Compute the signed distance on faces
    files = [ os.path.join(directories["aligned"],f) for f in os.listdir(directories["aligned"]) if "face.mesh" in f ]
    files = [ f for f in files if f.split("/")[-1][:3] + "_face.mesh" not in os.listdir(directories["signed"]) ]
    files.sort()
    run(signedDistance, files, maxi=4)

    # IV.2 - Morphing from template face to signedDistance faces
    files = [ os.path.join(directories["signed"],f) for f in os.listdir(directories["signed"]) if "face.mesh" in f ]
    files = [ [f, templates["morphing"]] for f in files if f.split("/")[-1][:3] + "_face.mesh" not in os.listdir(directories["morphed"]) ]
    files.sort()
    run(morph, files)


    ############################################################################
    # V - MASKS
    ############################################################################

    # V.1 - Generate the masks
    faces = [ os.path.join(directories["morphed"],f) for f in os.listdir(directories["morphed"]) if "face.mesh" in f]
    bones = [ os.path.join(directories["morphed"],f) for f in os.listdir(directories["morphed"]) if "bone.mesh"  in f]
    files = [[b,f] for b in bones for f in faces if f.split("/")[-1][:3]==b.split("/")[-1][:3]]
    files.sort()
    files = [f for f in files if f[0].split("/")[-1][:3] + "_mask.mesh" not in "".join(os.listdir(directories["masked"]))]
    run(generateMask, files)


    """
    f = "/Data/Facile3/raw/045_face.mesh"
    err1 = os.system( exe.boundingMesh + " -i " + f + " -o carved.mesh -r 51")
    err2 = os.system("blender --background --python blender_warp.py -- carved.mesh " + f + " warped.mesh")
    err3 = os.system(exe.mmgs + " -nr -nreg warped.mesh -o warped.o.mesh -hausd 0.002")
    newFaces = [ os.path.join(directories["remeshed"],f) for f in os.listdir(directories["remeshed"]) if "face.o.final.rewarped.mesh" in f]
    for f in newFaces:
        print f
    """



    """
    #MASSETERS AND MANDIBULES
    massFiles = [f for f in os.listdir(directories["raw"]) if "mass.mesh" in f]
    mandFiles = [f for f in os.listdir(directories["raw"]) if "mand.mesh" in f]
    massFiles.sort()
    mandFiles.sort()
    files = massFiles + mandFiles
    files = [f for f in files if not f[:-5] + ".R.mesh" in "".join(os.listdir(directories["muscles"]))]
    run(cut, files)

    align=False
    liste = os.listdir(directories["muscles"])
    files = [os.path.join(directories["muscles"],f) for f in liste if (".R." in f or ".L." in f)]
    files = [f for f in files  if f.split("/")[-1][:-5] + ".final." not in "".join(liste) and f.split("/")[-1][:-5] + ".warped." not in "".join(liste)]
    files.sort()

    if len(files)>0:
        if align:
            print str(len(files)) + " halved files to align"
            run(alignToTemplate,files)
            print str(len(files)) + " halved files to warp"
            run(warpWithBlender, files)

    files = [os.path.join(directories["muscles"],f) for f in os.listdir(directories["muscles"]) if ".warped" in f]
    run(signedDistance, files)

    #liste = os.listdir(directories["muscles"])
    #files = [os.path.join(directories["muscles"],f) for f in liste if (".R.final.mesh" in f or ".L.final.mesh" in f) and f[:-5] + ".final.o.mesh" not in "".join(liste)]
    #print str(len(files)) + " halved files to remesh"
    #if len(files)>0:
    #    run(remesh,files)
    """

    ############################################################################
    # RECONSTRUCTION
    ############################################################################
    """
    The new skull has to be reconstructed against the database:
    1 - Scale, align and warp it
    2 - Run a PCA on the distance to get a linear combination
    3 - Get the mean of the modified masks with a linear combination
    """
    recon = True

    if recon:

        os.chdir(directories["reconstruction"])

        unknownSkull = "/home/norgeot/unknownSkull.o2.mesh"#.mesh"
        skull = msh.Mesh(unknownSkull)

        #Scale
        if not os.path.exists("scaled.mesh"):
            scale = 0.0035
            skull.verts[:,:3] -= skull.center
            skull.verts[:,:3] *= scale
            skull.verts[:,:3] += [0.5,0.5,0.5]
            skull.write("scaled.mesh")

        #Align
        if not os.path.exists("aligned.mesh"):
            err = os.system(exe.align + " -i scaled.mesh " + templates["bone"] + " -d 0.1 -o 0.95 > /dev/null 2>&1")
            if err:
                raise FacileError("alignement failure")
            skull.applyMatrix(matFile = "mat_Super4PCS.txt")
            skull.applyMatrix(matFile = "mat_ICP.txt")
            skull.write("aligned.mesh")
            err = os.system(exe.pythonICP + " -s aligned.mesh -t " + templates["bone"] + " -m mat_pyICP.txt > /dev/null 2>&1")
            if err:
                pass
            else:
                skull.applyMatrix(matFile = "mat_pyICP.txt")
                skull.write("aligned.mesh")

        #Warp
        if not os.path.exists("warped.mesh"):
            os.system("cp " + templates["sphere"] + " ./sphere.mesh")
            err = os.system( exe.warping + " aligned.mesh -p -nit 150 -load 40 > /dev/null 2>&1" )
            if err:
                raise FacileError("Warping failure")
            warped = msh.Mesh("sphere.d.mesh")
            ext_ref = 2
            warped.tris = warped.tris[warped.tris[:,-1] != ext_ref]
            warped.tets = np.array([])
            warped.discardUnused()
            warped.write("warped.mesh")

        #Get the signed distance
        if not os.path.exists("signed.mesh"):
            signedDistance(os.path.join(directories["reconstruction"],"warped.mesh"))
            os.system("mv " + os.path.join(directories["signed"], "warped.mesh") + " ./signed.mesh")
            os.system("mv " + os.path.join(directories["signed"], "warped.sol")  + " ./signed.sol")

        #Morph the template skull to the signed distance
        if not os.path.exists("morphed.mesh"):
            morph([os.path.join(directories["reconstruction"], "signed.mesh"), templates["morphingSkull"]])
            os.system("mv " + os.path.join(directories["morphed"], "signed.mesh") + " ./morphed.mesh")
            os.system("mv " + os.path.join(directories["morphed"], "signed.sol")  + " ./morphed.sol")

        #Create the dirichlet files for every mask, from the mask to the unknown skull
        morphed = msh.Mesh("morphed.mesh")
        morphed.readSol()
        templateToUnknown = morphed.vectors
        files = [ os.path.join(directories["morphed"], f) for f in os.listdir(directories["morphed"]) if ".mesh" in f and "bone" in f]
        for f in files:
            skull = msh.Mesh(f)
            skull.readSol()
            templateToSkull = skull.vectors
            skullToUnknown = -templateToSkull + templateToUnknown




        #CHANGE FROM WRAPPING TO MORPHING
        sphere = msh.Mesh("sphere.d.mesh")

        #Get the closest warped surfaces
        files = [ f for f in os.listdir(directories["warped"]) if "bone" in f and "45" not in f ]
        files.sort()
        DISP = []
        for f in files:
            mesh = msh.Mesh(os.path.join(directories["warped"], f))
            disp = sphere.verts[:,:3] - mesh.verts[:,:3]
            DISP.append(disp)
        DISP = np.array(DISP)

        warped  = msh.Mesh("warped.mesh")
        unknown = sphere.verts[:,:3] - warped.verts[:,:3]

        #PCA
        n=25
        X = PCA(DISP, unknown, n)
        mesh=msh.Mesh()
        mesh.verts   = np.array([[x[0],x[1],x[2],0] for x in (sphere.verts[:,:3]-X)])
        mesh.scalars = np.array([ np.sqrt(np.sum([(x-y)**2 for x,y in zip(v1,v2)])) for v1,v2 in zip(unknown,X)])
        mesh.tris    = sphere.tris
        mesh.write("warped."+str(n)+".mesh")
        mesh.writeSol("warped."+str(n)+".sol")

        #Least squares on distances
        B = np.linalg.norm(unknown,axis=1)
        A = np.transpose(np.array([np.linalg.norm(D,axis=1) for D in DISP]))
        scalars = np.linalg.lstsq(A, B)
        coefs = np.array(scalars[0])
        print coefs
        proposed = np.sum([ c * d for c,d in zip(coefs, DISP) ], axis=0)
        mesh=msh.Mesh()
        mesh.verts   = np.array([[x[0],x[1],x[2],0] for x in (sphere.verts[:,:3]-proposed)])
        mesh.scalars = np.array([ np.sqrt(np.sum([(x-y)**2 for x,y in zip(v1,v2)])) for v1,v2 in zip(unknown,proposed)])
        mesh.tris    = sphere.tris
        mesh.write("warped.lstq.mesh")
        mesh.writeSol("warped.lstq.sol")

        #Reduced least squares with 5 approximations
        n = 5
        inds = np.abs(coefs).argsort()[-n:][::-1]
        A = np.transpose(np.array([np.linalg.norm(D,axis=1) for D in DISP[inds]]))
        scalars = np.linalg.lstsq(A, B)
        coefs = np.array(scalars[0])
        print coefs, inds
        proposed = np.sum([ c * d for c,d in zip(coefs, DISP[inds]) ], axis=0)
        mesh=msh.Mesh()
        mesh.verts   = np.array([[x[0],x[1],x[2],0] for x in (sphere.verts[:,:3]-proposed)])
        mesh.scalars = np.array([ np.sqrt(np.sum([(x-y)**2 for x,y in zip(v1,v2)])) for v1,v2 in zip(unknown,proposed)])
        mesh.tris    = sphere.tris
        mesh.write("warped.lstq5.mesh")
        mesh.writeSol("warped.lstq5.sol")




        #Get the mean warped bone
        """
        meanSkull = msh.Mesh()
        meanSkull.verts = np.zeros(shape = (len(warped.verts), 4))
        meanSkull.verts[:,:3] = sphere.verts[:,:3] - np.mean(DISP, axis=0)
        meanSkull.tris = warped.tris
        meanSkull.tris[:,-1]=1
        meanSkull.write("/home/norgeot/templateMorphingSkull.mesh")
        """


        #Morph the unknown skull bone to the templateMesh
        #if not os.path.exists(os.path.join(directories["signed"],"templateMorphingSkull.mesh")):


        """
        #Tests pour l'élasticité
        files = [ os.path.join(directories["masked"], f) for f in os.listdir(directories["masked"]) if ".o.mesh" in f ]
        for f in files:
            mesh2 = msh.Mesh(f)
            n = len(mesh2.verts)
            mesh2.tris = mesh2.tris[mesh2.tris[:,-1]==1]
            mesh2.tets = np.array([])
            mesh2.discardUnused()

            mesh2.vectors = np.zeros((n,3))
            mesh2.vectors[:len(mesh.verts)] = mesh2.verts[:,:3] - warped.verts[:,:3] # A remplacer par le resultat du morphing
            print len(mesh2.vectors)
            mesh2.writeSol( os.path.join(directories["masked"], f.split("/")[-1].split(".")[0] + ".o.sol") )
        """

        #On garde les coefficients du least squares pour faire derriere la moyenne des morphing obtenus
