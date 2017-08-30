import numpy as np
import time
import matplotlib.pyplot as plt

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
    plt.figure()
    plt.plot(np.log10(eVal))
    plt.grid()
    plt.show()
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
