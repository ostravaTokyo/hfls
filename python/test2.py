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


K_modif_0_v1 = load_matrix_basic("../data/dump_K_modif_v1_0.txt",False,True, 1)
K_modif_1_v1 = load_matrix_basic("../data/dump_K_modif_v1_1.txt",False,True, 1)
                                                                      
K_modif_0_v2 = load_matrix_basic("../data/dump_K_modif_v2_0.txt",False,True, 1)
K_modif_1_v2 = load_matrix_basic("../data/dump_K_modif_v2_1.txt",False,True, 1)
                                                                      
K_modif_0_v3 = load_matrix_basic("../data/dump_K_modif_v3_0.txt",False,True, 1)
K_modif_1_v3 = load_matrix_basic("../data/dump_K_modif_v3_1.txt",False,True, 1)


K_rs_0 = load_matrix_basic("../data/dump_K_rs_0.txt",False,False, 1)
K_rs_1 = load_matrix_basic("../data/dump_K_rs_1.txt",False,False, 1)

K_rr_0 = load_matrix_basic("../data/dump_K_rr_0.txt",False,True, 1)
K_rr_1 = load_matrix_basic("../data/dump_K_rr_1.txt",False,True, 1)


K_ss_0 = load_matrix_basic("../data/dump_K_ss_0.txt",False,True, 1)
K_ss_1 = load_matrix_basic("../data/dump_K_ss_1.txt",False,True, 1)


InvKrrKrs_0 = load_matrix_basic("../data/dump_InvKrrKrs_0.txt",False,False, 1)
InvKrrKrs_1 = load_matrix_basic("../data/dump_InvKrrKrs_1.txt",False,False, 1)

InvKrrKrs_0_p = np.linalg.solve(K_rr_0,K_rs_0)
InvKrrKrs_1_p = np.linalg.solve(K_rr_1,K_rs_1)

#dump_Precond_0.txt
P_0 = load_matrix_basic("../data/dump_Precond_0.txt",False,False, 1)
P_1 = load_matrix_basic("../data/dump_Precond_1.txt",False,False, 1)


#dump_KsrInvKrrKrs_0
KsrInvKrrKrs_0 = load_matrix_basic("../data/dump_KsrInvKrrKrs_0.txt",False,False, 1)
KsrInvKrrKrs_1 = load_matrix_basic("../data/dump_KsrInvKrrKrs_1.txt",False,False, 1)


KsrInvKrrKrs_0_p = np.dot(K_rs_0.T,np.linalg.solve(K_rr_0,K_rs_0))
KsrInvKrrKrs_1_p = np.dot(K_rs_0.T,np.linalg.solve(K_rr_1,K_rs_1))

#plt.spy(K_modif_0_v3);plt.show()
