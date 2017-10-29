
import numpy as np
#from scipy.linalg import solve
from scipy import sparse
#from scipy import linalg
import myModul as mM
import config_espreso_python

#import pylab as plt  


n_clus          = 1
n_subPerClust   = 2

path0 = '/data_space/WorkSpace/htfeti_app/gitRepo/hfls/data_1/'

mat_K       = []
#mat_Krs       = []
#mat_Krr       = []
#mat_Kss       = []
mat_S       = []
mat_Kreg    = []
mat_B0      = []
mat_B1      = []
mat_R       = []
vec_f       = []
vec_weight  = np.loadtxt(path0+"/dump_weigth.txt")
for i in range(n_clus): 
    mat_K.append([])
#    mat_Krr.append([])
#    mat_Kss.append([])
#    mat_Krs.append([])
    mat_S.append([])
    mat_Kreg.append([])
    mat_B0.append([])
    mat_B1.append([])
    mat_R.append([])
    vec_f.append([])
#    vec_weight.append([])
    for j in range(n_subPerClust):  
        mat_K[i].append(mM.load_matrix0(path0,"dump_K_","",str(j),True,True,1))
#        mat_Krr[i].append(mM.load_matrix0(path0,"dump_Krr_","",str(j),True,True,1))
#        mat_Krs[i].append(mM.load_matrix0(path0,"dump_Krs_","",str(j),True,False,1))
#        mat_Kss[i].append(mM.load_matrix0(path0,"dump_Kss_","",str(j),True,True,1))
        mat_S[i].append(mM.load_matrix0(path0,"dump_S_","",str(j),True,False,1))
        mat_Kreg[i].append(mM.load_matrix0(path0,"dump_K_reg_","",str(j),True,True,1))
        mat_B0[i].append(mM.load_matrix0(path0,"dump_Bc_","",str(j),True,False,1))
        mat_R[i].append(mM.load_matrix0(path0,"dump_R_","",str(j),True,False,1))
        vec_f[i].append(mM.load_matrix0(path0,"dump_rhs_","",str(j),False,False,1))
        mat_B1[i].append(mM.load_matrix0(path0,"dump_Bf_","",str(j),True,False,1))
        #vec_weight[i].append(mM.load_vector(path,'weight',i,j))
        
    
for i in range(n_clus):
    for j in range(n_subPerClust):
        if (i==0 and j==0):
            B0      = mat_B0[0][0].copy()  
            print('1')  
            B1      = mat_B1[0][0].copy()      
            print('2')  
            K       = mat_K[0][0]   
            print('3')  
            Kreg    = mat_Kreg[0][0]
            print('4')  
            R       = mat_R[0][0]
            print('5')  
            f       = vec_f[0][0]
            print('6')  
            #weight  = vec_weight[0][0]
#            diagR   = np.sum(mat_R[0][0]*mat_R[0][0],axis=1)
            print('7')  
            S       = mat_S[0][0]

        else:  
            B0      = sparse.hstack((B0,mat_B0[i][j]))
            print('-1')  
            B1      = sparse.hstack((B1,mat_B1[i][j]))
            print('-2')  
            K       = sparse.block_diag((K,mat_K[i][j]))
            print('-3')  
            Kreg    = sparse.block_diag((Kreg,mat_Kreg[i][j]))
            print('-4')  
            R       = sparse.block_diag((R,mat_R[i][j]))            
            print('-5')  
            f       = np.concatenate((f,vec_f[i][j]))        
            print('-6')  
            #weight  = np.concatenate((weight,vec_weight[i][j]))
#            diagR   = np.concatenate((diagR,np.sum(mat_R[i][j]*mat_R[i][j],axis=1)))           
            print('-7')  
            S       = sparse.block_diag((S,mat_S[i][j]))
            
        
#diagR   = diagR
B0      = B0.tocsc()
if n_clus*n_subPerClust==1:
    R       = sparse.csc_matrix(R)
B1      = B1.tocsc()    
K       = K.tocsc()     
Kreg    = Kreg.tocsc()  
R       = R.tocsc()     
weight = 1
            
#mat_K       = []
#mat_Kreg    = []
#mat_B0      = []
#mat_B1      = []
#mat_R       = []
#vec_f       = []
#vec_weight  = []            
###############################################################################
####################### FETI PREPROCESSING ####################################
###############################################################################            
#np.hstack;np.vstack;np.column_stack;np.row_stack


#FETI=False

#if FETI:
    #B = sparse.vstack((B0 ,B1 ))
    #u,lam = mM.feti(K,Kreg,f,B1,R,weight)
#else:    
    #u,lam = mM.hfeti(K,Kreg,f,B0,B1,R,weight)

#B = sparse.vstack((B0 ,B1 ))
#u,lam = mM.feti(K,Kreg,f,B,R,weight)

conf = config_espreso_python


#u,lam = mM.feti(K,Kreg,f,B1,R,diagR,weight)
uHDP,lamH = mM.hfeti(K,Kreg,f,B0,B1,R,S, vec_weight)

#print("u")
#for i in range(uHDP.shape[0]):
#    print("%d %e" % (i+1,uHDP[i]))


#print("u\n", uHDP)
#print("lamH\n",lamH)

#conf.iterative_Kplus=False
#uHDP,lamH = mM.hfeti(K,Kreg,f,B0,B1,R,weight)
#
#methodToImproveSolByIterMethod      = 'cg_dx'
#conf.precondPrimalSystem            = 'diag' 
#conf.iterative_Kplus=True
#uHSP,lamH = mM.hfeti(K,Kreg,f,B0,B1,R,weight)



#ndu = np.linalg.norm(uHSP-uHDP)
#nu = np.linalg.norm(uHDP)
#print('||uHDP-uHSP||/||uHDP||       = ',ndu/nu)

#ndu = np.linalg.norm(u-uH)
#nu = np.linalg.norm(u)
#ndlam = np.linalg.norm(lam-lamH)
#nlam = np.linalg.norm(lam)
#print('||u-uH||/||u||       = ',ndu/nu)
#print('||lam-lamH||/||lam|| = ',ndlam/nlam)
