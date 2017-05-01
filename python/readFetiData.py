
import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import pylab as plt
#
#
def load_matrix(path,str0,i,j,makeSparse,makeSymmetric,offset): 
    pathToFile = path+'/'+str(i)+'/'+str0+str(j)+'.txt' 
    tmp = np.loadtxt(pathToFile, ndmin=2)
    if (tmp.shape[0]==1):
        tmp = []
    else:
        n = np.int32(tmp[0,0])   
        m = np.int32(tmp[0,1])
        I = tmp[1::,0]-offset;
        J = tmp[1::,1]-offset;
        V = tmp[1::,2]
#    
        print(str0,i,j)
        if (makeSymmetric):
            logInd = J>I; 
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
#
    return tmp

path0 = "../data"


Bc_List = []
Bct_list = []

K_newList = []
R_newList = []

Gc_newList = []
Gc_List = []

K_List = []

K_reg_List = []
Fc_List = []
BcKplus_List = []


nSub = 8




for i in range(nSub):

    K_new  = load_matrix(path0,"dump_K_","",str(i),True,True,0)
    K_orig = load_matrix(path0,"K","",str(i),True,False,1) 
    del_ = np.linalg.norm(K_new.toarray() - K_orig.toarray())
    n_k_orig = np.linalg.norm(K_orig.toarray()) 

    K_List.append(K_new.toarray())

    K_reg = load_matrix(path0,"dump_K_reg_","",str(i),False,True,0)
    K_reg_List.append(K_reg)


    Fc = load_matrix(path0,"dump_Fc_","",str(i),False,False,0)
    Fc_List.append(Fc)

    K_newList.append(K_new)

    print "|K_new - K_orig| / | K_orig | :    %e, | K_orig | = %e " % (del_, n_k_orig) 

    R_new  = load_matrix(path0,"dump_R_","",str(i),True,False,0)
    R_newList.append(R_new)
    R_orig = load_matrix(path0,"R1","",str(i),True,False,1) 
    del_ = np.linalg.norm(R_new.toarray() - R_orig.toarray())
    n_k_orig = np.linalg.norm(R_orig.toarray()) 
    print "|R_new - R_orig| / | R_orig | :    %e, | R_orig | = %e " % (del_, n_k_orig) 


    
    Bc_new  = load_matrix(path0,"dump_Bc_","",str(i),True,False,0)
    Bc_List.append(Bc_new)
    Bc_orig = load_matrix(path0,"B0","",str(i),True,False,1) 
    del_ = np.linalg.norm(Bc_new.toarray() - Bc_orig.toarray())
    n_k_orig = np.linalg.norm(Bc_orig.toarray()) 
    print "|Bc_new - Bc_orig| / | Bc_orig | : %e, | Bc_orig | = %e " % (del_, n_k_orig) 
    Bct_new  = load_matrix(path0,"dump_Bc_dense_","",str(i),True,False,0)
    Bct_list.append(Bct_new)

    Bc = Bc_List[0].toarray();
    Rc = R_newList[0].toarray();
    
    Gc_List.append(np.dot(Bc, Rc))
    Gc_newList.append(load_matrix(path0,"dump_Gc_","",0,False,False,0))


    BcKplus  = load_matrix(path0,"dump_BcKplus_","",str(i),False,False,0)
    BcKplus_List.append(BcKplus)
#    
#    
#    Bf_new  = load_matrix(path0,"dump_Bf_","",str(i),True,False)
#    Bf_orig = load_matrix(path0,"B1","",str(i),True,False) 
#    del_ = np.linalg.norm(Bf_new.toarray() - Bf_orig.toarray())
#    n_k_orig = np.linalg.norm(Bf_orig.toarray()) 
#    print "\t\t|Bf_new - Bf_orig| / | Bf_orig | : \n%e,    | Bf_orig | = %e " % (del_, n_k_orig) 

#K_orig = load_matrix(path0,"K","",str(i),True,False) 
#R = load_matrix(path0+'1',"dump_Kreg_","",str(i),True,True) 

Fc_python_List = []
for i in range(nSub):
    Bc = Bct_list[i].toarray()
    indBc = np.abs(Bc).sum(axis=1)>0
    Bc = Bc[indBc,:]
    BcKplus = BcKplus_List[i]
    indBcKplus = np.abs(BcKplus).sum(axis=1)>0
    BcKplus = BcKplus[indBcKplus,:] 
    BcKplus_python = np.linalg.solve(K_reg_List[i],Bc.T)
    Fc_i = np.dot(Bc,BcKplus_python)
    Fc_python_List.append(Fc_i)

for ii in range(nSub):
 
    ttt = Gc_List[ii][np.abs(Gc_List[ii]).sum(axis=1)>0,:] - Gc_newList[ii]
    print( np.linalg.norm(ttt))


for ii in range(nSub):
    ddd0 = np.linalg.norm(Fc_python_List[ii] - Fc_List[ii])
    ddd1 = np.linalg.norm(Fc_python_List[ii])
    print( ddd0 / ddd1 )
