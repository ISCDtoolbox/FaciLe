import numpy as np
import os
import cPickle as pickle
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../tools"))
import msh

if __name__=="__main__":
    """
    Runs a PCA on all the input files
    Usage:
        -- python PCA.py -PCA *.sol : Runs the PCA
        -- python PCA.py -REC object.sol n: Reconstruct object.sol from the PCA with n coefficients
    """

    if sys.argv[1] == "-PCA":
        ###################################################################
        # Running the PCA
        ###################################################################
        print "Running the PCA on input files"

        # Files reading
        print "- Reading the vector fields"
        files = [f for f in sys.argv[2:] if ".sol" in f]
        files.sort()
        data = []
        for f in files:
            mesh = msh.Mesh()
            mesh.readSol(f)
            data.append(mesh.vectors)
        data = np.array(data)
        data = np.array([d - np.mean(d,axis=0) for d in data])

        # Scalar products
        print "- Computing the scalar products"
        A = np.array([[np.sum(np.multiply(x,y)) for x in data] for y in data])

        # Eigen values and vectors
        print "- Computing the eigen values and vectors"
        eVal, eVec = np.linalg.eig(A)
        eigenVectors = np.array([vec for (val,vec) in sorted(zip(eVal,eVec))[::-1]])
        eigenValues  = np.array([val for (val,vec) in sorted(zip(eVal,eVec))[::-1]])

        # Getting the PCA coefficients
        B = [ np.sum([eigenVectors[i,j] * data[i] for i in range(len(data))],axis=0) for j in range(len(data))]
        BB = [b/np.sqrt(np.sum(np.multiply(b,b))) for b in B]
        print np.array(BB).shape

        # Saving
        print "- Dumping the data"
        print data.shape
        with open("PCA_data", 'wb') as f:
            pickle.dump(data, f)
            pickle.dump(BB, f)

    elif sys.argv[1] == "-REC":
        ###################################################################
        # Reconstructing from a previous PCA
        ###################################################################
        print "Reconstructing the given object from the PCA"

        # Opening the file
        print "- Opening the file"
        mesh = msh.Mesh()
        mesh.readSol(sys.argv[2])
        unknown = mesh.vectors - np.mean(mesh.vectors, axis=0)

        # Getting the PCA data
        print "- Loading the data"
        with open("PCA_data", 'rb') as f:
            data = pickle.load(f)
            BB   = pickle.load(f)

        # Computing the reconstruction coefficients
        print "- Computing the coefficients"
        sc = float(0.5/len(data))
        alpha = np.array([[sc*np.sum(np.multiply(d,bb)) / np.sum(np.multiply(bb,bb)) for d in np.append(data,[unknown],axis=0)] for bb in BB])

        #Combinaison lineaire
        """
        print "Combinaison lineaire:"
        L = np.zeros((nData + 1, nPts, 3))
        for i in range(nData):
            print '{0}\r'.format("  " + str(i+1) + "/" + str(nData)),
            for j in range(nData + 1):
                L[j] += alpha[i,j] * BB[i]
        print "\n"
        #L = [ np.sum([a*bb for a in alpha[i]],axis=0) for i,bb in enumerate(BB)]
        """

        #N = [np.sum([alpha[i,j] * BB[i] for i in range(nData-3,nData)],axis=0) for j in range(nData + 1)]
        #print(N[0])
        print "- Linear combination"
        N = [np.sum([a*BB[i] for a in alpha[i]], axis=0) for i in range(int(sys.argv[3]))]
        X = np.sum(N,axis=0)
        mesh=msh.Mesh()
        mesh.vectors = X
        mesh.writeSol("toto.sol")
