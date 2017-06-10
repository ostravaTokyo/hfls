import numpy as np 
from scipy import sparse
import scipy.sparse.linalg as spla
import pylab as plt
#
#

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
nSub = 8




if 1:

    K = []
#    K_reg = []
    Fc = []
    R = []
    Bc = []
    BcT_dense = []
    Gc = []

    Gc = []
    Fc = []


    Bc_nonzRow = []
                 
    BcKplus = []
    BcKplus = []
    BcKplus_tmp = []

#    BcK_dense = []

    K_UT = []

#    Lumped = []
#    Lumped = []




    for i in range(nSub):

        K.append(load_matrix(path0,"dump_K_","",str(i),False,True,1)) 
        K_UT.append(load_matrix(path0,"dump_K_","",str(i),False,False,1)) 
#        K_reg.append(load_matrix(path0,"dump_K_reg_","",str(i),False,True,1)) 
        Fc.append(load_matrix(path0,"dump_Fc_","",str(i),False,False,1))
        R.append(load_matrix(path0,"dump_R_","",str(i),False,False,1)) 
        Bc.append(load_matrix(path0,"dump_Bc_","",str(i),False,False,1))
#        Lumped.append(load_matrix(path0,"dump_Lumped_","",str(i),False,False,1))
        BcT_dense.append(load_matrix(path0,"dump_BcT_dense_","",str(i),False,False,1))
        Gc.append(load_matrix(path0,"dump_Gc_","",str(i),False,False,1))

        indBc = np.abs(Bc[i]).sum(axis=1)>0
        Bc_nonzRow.append( Bc[i][indBc,:])
#        Fc.append( np.dot(Bc_nonzRow[i], np.linalg.solve(K_reg[i],Bc_nonzRow[i].T)))
#        Lumped.append( np.dot(Bc_nonzRow[i], np.dot(K[i],Bc_nonzRow[i].T)))



    #    BcKplus = BcKplus_List[i]

#        BcK_dense.append(load_matrix(path0,"dump_BcK_dense_","",str(i),False,False,1))
#        BcK_dense.append(np.dot(K[i],Bc_nonzRow[i].T).T)
        
        Gc.append(np.dot(Bc[i], R[i]))

#        BcKplus.append(load_matrix(path0,"dump_BcKplus_","",str(i),False,False,1))
#        BcKplus.append(np.linalg.solve(K_reg[i],Bc_nonzRow[i].T).T)
#        BcKplus_tmp.append(np.linalg.solve(K_reg[i],Bc[i].T).T)

#        iK_K = np.linalg.solve(K_reg[i],K[i])
#        K_iK_K = np.dot(K[i],iK_K)
#        del_ = np.linalg.norm(K_iK_K - K[i]  ) / np.linalg.norm(K[i])
#        print(del_)
#
        print(' ...%d '%(i))
    Gc_clust = load_matrix(path0,"dump_Gc_clust_","",str(0),False,False,1)
    Ac_clust = load_matrix(path0,"dump_Ac_clust_","",str(0),False,True,1)
    Fc_clust = load_matrix(path0,"dump_Fc_clust_","",str(0),False,True,1)

    Ac_clust_python = np.hstack((Fc_clust,Gc_clust))

    Z = np.zeros((Gc_clust.shape[1],Ac_clust_python.shape[1]))
    print ( Z.shape)
    Ac_clust_python = np.vstack((Ac_clust_python,Z))


for i in range(nSub - 1):
    if (i == 0):
        Bc_g = np.hstack((Bc[0],Bc[1]))
    else:
        Bc_g = np.hstack((Bc_g,Bc[i+1]))



#Bc_g = np.hstack((Bc_g,Bc[2]))
#Bc_g = np.hstack((Bc_g,Bc[2])) 
#BcT_dense = load_matrix(path0,"dump_BcT_dense_","",str(0),True,True,1) 
Fc_clust = load_matrix(path0,"dump_Fc_clust_","",str(0),True,True,1)
Ac_clust = load_matrix(path0,"dump_Ac_clust_","",str(0),True,True,1)
GcTGc = load_matrix(path0,"dump_GcTGc_clust_","",str(0),False,True,1) 
ker_GcTGc = load_matrix(path0,"dump_kerGc_","",str(0),False,False,1) 
ker_Ac = load_matrix(path0,"dump_ker_Ac_","",str(0),False,False,1) 
#KpBcT0 = load_matrix(path0,"dump_KplusBcT_","",str(0),False,False,1) 
#KpBcT1 = load_matrix(path0,"dump_KplusBcT_","",str(1),False,False,1) 


dFc_eig = load_matrix(path0,"dump_Fc_clust_","",str(444),False,False,1)
#dFc_svd = load_matrix(path0,"dump_Fc_clust_","",str(555),False,False,1)
dAc_eig = load_matrix(path0,"dump_Ac_clust_","",str(444),False,False,1)
#dAc_svd = load_matrix(path0,"dump_Ac_clust_","",str(555),False,False,1)



plt.subplot(1,3,1)
if GcTGc.shape[0] < 100:
    markersize_ = 3
else:
    markersize_ = 0.7 
plt.spy(GcTGc, markersize=markersize_)
plt.xlabel("nnz = %d" % (GcTGc.nonzero()[0].shape[0]))
plt.subplot(1,3,2)
if Fc_clust.shape[0] < 100:
    markersize = 3
else:
    markersize = 0.7 
plt.spy(Fc_clust, markersize=markersize_)
plt.xlabel("nnz = %d" % (Fc_clust.nonzero()[0].shape[0]))
plt.subplot(1,3,3)
if Ac_clust.shape[0] < 100:
    markersize_ = 3
else:          
    markersize_ = 0.7 
plt.spy(Ac_clust, markersize=markersize_) 
plt.xlabel("nnz = %d" % (Ac_clust.nonzero()[0].shape[0]))
plt.show()

#Bc_from_Rt = []
#for i in range(1,14):
#    Bc_from_Rt.append( load_matrix(path0,"dump_Bc_from_Rt_","",str(i),False,False,1) )
#

# Gc_ = load_matrix(path0,"dump_Gc_i_","",str(0),False,False,1)




#BcT_dense = load_matrix(path0,"dump_BcT_dense_","",str(0),True,True,1) 


#K_test= []
#Kplus_K_test = []
#K_Kplus_K_test = []
#K_reg_test = []
#for i in range(2): 
#    K_test.append(load_matrix(path0,"dump_K_dense_","",str(i),False,False,1))
#    K_reg_test.append(load_matrix(path0,"dump_K_reg_","",str(i),False,True,1)) 
#    Kplus_K_test.append(load_matrix(path0,"dump_Kplus_K_","",str(i),False,False,1))
#    K_Kplus_K_test.append(load_matrix(path0,"dump_K_Kplus_K_","",str(i),False,False,1))
#

































#plt.spy(Fc_clust,markersize = .8);plt.show()

#Gc_ = load_matrix(path0,"dump_Gc_i_","",str(0),True,True,1)



#r = sparse.csgraph.reverse_cuthill_mckee(Ac_clust.tocsr(), symmetric_mode=True)
#Ac_clust = Ac_clust.toarray() 
##
#P,L,U= scipy.linalg.lu(Ac_clust)
#nnz0 = L.nonzero()[0].shape[0] +  U.nonzero()[0].shape[0]
#
#
##
##
#AcR = Ac_clust[np.ix_(r,r)]
#PR,LR,UR = scipy.linalg.lu(AcR)
#nnzR = LR.nonzero()[0].shape[0] +  UR.nonzero()[0].shape[0]
##
##
#plt.subplot(2,2,1)
#plt.spy(L,markersize=0.1);
#plt.subplot(2,2,2)
#plt.spy(U,markersize=0.1);
#plt.subplot(2,2,3)
#plt.spy(LR,markersize=0.1);
#plt.subplot(2,2,4)
#plt.spy(UR,markersize=0.1);

#print ("nnz = %d, nnz(reordered) = %d ") % (nnz0, nnzR) 


#plt.show()

#ker_Ac = load_matrix(path0,"dump_ker_Ac_","",str(0),False,True,1) 
#ker_GcTGc = load_matrix(path0,"dump_ker_GcTGc_","",str(0),False,True,1) 
#R0 = load_matrix(path0,"dump_R_","",str(0),False,True,1) 

#Gc_H = np.dot(GcTGc.toarray(),ker_GcTGc)

#r = sparse.csgraph.reverse_cuthill_mckee(Ac_clust.tocsr(), symmetric_mode=True)
#Ac = Ac_clust.toarray()[np.ix_(r,r)]
#plt.subplot(1,2,1)
#plt.spy(Ac_clust ,markersize = 2.0)
#plt.subplot(1,2,2)
#plt.spy(Ac,markersize = 0.125)



#Fc_python_List = []

#if 0:
#    Fc_clust_python = np.zeros((Bct_list[i].shape[0], Bct_list[i].shape[0]))
#    for i in range(nSub):
#        Bc = Bct_list[i].toarray()
#        indBc = np.abs(Bc).sum(axis=1)>0
#        Bc_red = Bc[indBc,:]
#        BcKplus = BcKplus_List[i]
#
#        Bf = Bf_List[i].toarray()
#        indBf = np.abs(Bf).sum(axis=1)>0
#        Bf_red = Bf[indBf,:]
#
#        Rc = RList[i].toarray()
#        
#        
#
#        if (i == 0):
#            Gf_clust_python = np.dot(Bf,Rc)
#            Gc_clust_python = np.dot(Bc,Rc)
#        else:
#            Gf_clust_python = np.hstack((Gf_clust_python,np.dot(Bf,Rc)))
#            Gc_clust_python = np.hstack((Gc_clust_python,np.dot(Bc,Rc)))
#        indBcKplus = np.abs(BcKplus).sum(axis=1)>0
#        BcKplus = BcKplus[indBcKplus,:] 
#        BcKplus_python = np.linalg.solve(K_reg_List[i],Bc_red.T)
#        BcKplus_ = np.linalg.solve(K_reg_List[i],Bc.T)
#        Fc_i = np.dot(Bc_red,BcKplus_python)
#        Fc_clust_python += np.dot(Bc,BcKplus_)
#        Fc_python_List.append(Fc_i)
#
#    for ii in range(nSub):
#     
#        ttt = Gc_List[ii][np.abs(Gc_List[ii]).sum(axis=1)>0,:] - GcList[ii]
#        print np.linalg.norm(ttt)
#
#
#    for ii in range(nSub):
#        ddd0 = np.linalg.norm(Fc_python_List[ii] - Fc_List[ii])
#        ddd1 = np.linalg.norm(Fc_python_List[ii])
#        print "|Fc_python - Fc_myAp|/|Fc_python|",ddd0 / ddd1  
#
#
#    Fc_clust =  load_matrix(path0,"dump_Fc_clust_","",0,False,True,1)
#    Gc_clust =  load_matrix(path0,"dump_Gc_clust_","",0,False,False,1)
#    Gf_clust =  load_matrix(path0,"dump_Gf_clust_","",0,False,False,1)
#    Ac_clust =  load_matrix(path0,"dump_Ac_clust_","",0,False,True,1)
#    Ac_clust_python = np.hstack((Fc_clust_python,Gc_clust_python))
#
#    Z = np.zeros((Gc_clust_python.shape[1],Ac_clust.shape[1]))
#    print ( Z.shape)
#    Ac_clust_python = np.vstack((Ac_clust_python,Z))
#
#
#    ddd0 = np.linalg.norm(Fc_clust - Fc_clust_python)
#    ddd1 = np.linalg.norm(Fc_clust)
#    print "|Fc_clust_python - Fc_clust_myAp|/|Fc_clust_python|",ddd0 / ddd1  
#
#    ddd0 = np.linalg.norm(Gc_clust - Gc_clust_python)
#    ddd1 = np.linalg.norm(Gc_clust)
#    print "|Gc_clust_python - Gc_clust_myAp|/|Gc_clust_python|",ddd0 / ddd1  
#    
#    ddd0 = np.linalg.norm(Gf_clust - Gf_clust_python)
#    ddd1 = np.linalg.norm(Gf_clust)
#    print "|Gf_clust_python - Gf_clust_myAp|/|Gf_clust_python|",ddd0 / ddd1  
#
#    ddd0 = np.linalg.norm(Ac_clust - Ac_clust_python)
#    ddd1 = np.linalg.norm(Ac_clust)
#    print "|Ac_clust_python - Ac_clust_myAp|/|Ac_clust_python|",ddd0 / ddd1  
#
#K = []



#plt.subplot(1,2,1)
#plt.spy(Gf_clust_python,markersize=1)
#plt.subplot(1,2,2)
#plt.spy(Gf_clust,markersize=1)
#plt.show()
