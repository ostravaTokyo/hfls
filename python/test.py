import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import pylab as plt
from scipy.linalg import block_diag
#
#

nSub = 2
def load_matrix_basic(pathToFile,makeSparse,makeSymmetric, offset):
    f0 = open(pathToFile).readlines()
    firstLine = f0.pop(0) #removes the first line
    tmp = np.zeros((len(f0),3), dtype = float)
    for i in range(len(f0)):
        line = f0[i]
        k = line.split()
        tmp[i,0] = float(k[0])
        tmp[i,1] = float(k[1])
        tmp[i,2] = float(k[2])


    if (tmp.shape[0]==1):
        tmp = []
    else:
        n = np.int32(tmp[0,0])   
        m = np.int32(tmp[0,1])
        I = tmp[1::,0]-offset;
        J = tmp[1::,1]-offset;
        V = tmp[1::,2]
#    
#        print str0,i,j
        if (makeSymmetric):
            logInd = J != I; 
            I = np.concatenate((I,J[logInd]))
            J = np.concatenate((J,I[logInd]))
            V = np.concatenate((V,V[logInd]))    
        
        if (makeSparse):
            tmp = sparse.csc_matrix((V,(I,J)),shape=(n,m)).tocoo()
        else:
            if (m==1):
                tmp = V
            else:                
                tmp = sparse.csc_matrix((V,(I,J)),shape=(n,m)).toarray()
    return tmp

def load_matrix(path,str0,i,j,makeSparse,makeSymmetric,offset): 
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' #
    tmp = load_matrix_basic(pathToFile,makeSparse,makeSymmetric,offset) 
    return tmp

path0 = "../data"


