import numpy as np
import os
import cPickle as pickle
import matplotlib.pyplot as plt

def read_data(files):
    data = []
    for filename in files:
        values = []
        with open(filename,"r") as f:
            for l in f.readlines()[8:-2]:
                values.append([float(x) for x in l.split()])
        data.append(np.array(values))
    data = np.array([d-np.mean(d) for d in data])
    return np.array(data)

def plot_matrix(mat, file):
    plt.figure()
    plt.matshow(mat)
    plt.colorbar()
    plt.savefig(file)

if __name__ == "__main__":
    #Lecture des fichiers
    print "Lecture des donnees:"
    directory = "data/"
    files     = [os.path.join(directory,f) for f in os.listdir(directory) if ".sol" in f]
    files.sort()
    data = read_data(files)
    unknown = data[10]
    a = list(data[:10])
    a.extend(list(data[11:]))
    data = np.array(a)
    nData = len(data)
    nPts  = len(data[0])
    print "  OK\n"

    #Matrice des produits scalaires
    A = np.array([[np.sum(np.multiply(x,y)) for x in data] for y in data])
    plot_matrix(A,"A.png")

    eVal, eVec = np.linalg.eig(A)
    eigenVectors = np.array([vec for (val,vec) in sorted(zip(eVal,eVec))[::-1]])
    eigenValues  = np.array([val for (val,vec) in sorted(zip(eVal,eVec))[::-1]])

    B = [ np.sum([eigenVectors[i,j] * data[i] for i in range(nData)],axis=0) for j in range(nData)]

    BB = [b/np.sqrt(np.sum(np.multiply(b,b))) for b in B]

    sc = float(0.5/nData)
    print sc
    alpha = np.array([[sc*np.sum(np.multiply(d,bb)) / np.sum(np.multiply(bb,bb)) for d in np.append(data,[unknown],axis=0)] for bb in BB])
    plot_matrix(alpha,"alpha.png")

    #Combinaison lineaire

    print "Combinaison lineaire:"
    L = np.zeros((nData + 1, nPts, 3))
    for i in range(nData):
        print '{0}\r'.format("  " + str(i+1) + "/" + str(nData)),
        for j in range(nData + 1):
            L[j] += alpha[i,j] * BB[i]
    print "\n"
    #L = [ np.sum([a*bb for a in alpha[i]],axis=0) for i,bb in enumerate(BB)]

    #N = [np.sum([alpha[i,j] * BB[i] for i in range(nData-3,nData)],axis=0) for j in range(nData + 1)]
    #print(N[0])
    N = [np.sum([a*BB[i] for a in alpha[i]], axis=0) for i in range(10)]
    #print(N[0])

    with open("template3D.sol","w") as f:
        f.write("MeshVersionFormatted 1\n")
        f.write("Dimension 3\n")
        f.write("SolAtVertices\n")
        f.write(str(len(B[0])) + "\n")
        f.write("1 2\n")
        X = np.sum(N,axis=0)
        #X = L[5]
        for t in X:
            f.write(" ".join([str(x) for x in t]) + "\n")
