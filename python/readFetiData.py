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





if 1:

    K = []
    K_reg = []
    Fc = []
    R = []
    Rf = []
    Bc = []
    Bf = []
    BcT_dense = []
    Gc = []
    Gf = []
    Gf_p = []

    Gc = []
    Fc_p = []

    rhs = []
    xx = []
    Kplus_f_test = []
    KplusBcT_p = []


    Bc_nonzRow = []
                 
    KplusBcT = []
    BcKplus_tmp = []

#    BcK_dense = []

    K_UT = []

#    x_out = []
#    x_out_p = []
#    Lumped = []
#    Lumped = []




    for i in range(nSub):

        K.append(load_matrix(path0,"dump_K_","",str(i),False,True,1)) 
        K_UT.append(load_matrix(path0,"dump_K_","",str(i),False,False,1)) 
        K_reg.append(load_matrix(path0,"dump_K_reg_","",str(i),False,True,1)) 
        Fc.append(load_matrix(path0,"dump_Fc_","",str(i),False,False,1))
        R.append(load_matrix(path0,"dump_R_","",str(i),False,False,1)) 
        Rf.append(load_matrix(path0,"dump_Rf_","",str(i),False,False,1)) 
        Bc.append(load_matrix(path0,"dump_Bc_","",str(i),False,False,1))
        Bf.append(load_matrix(path0,"dump_Bf_","",str(i),False,False,1))
        Gf_p.append(np.dot(Bf[i],Rf[i]))
#        Lumped.append(load_matrix(path0,"dump_Lumped_","",str(i),False,False,1))
        BcT_dense.append(load_matrix(path0,"dump_BcT_dense_","",str(i),False,False,1))
        Gc.append(load_matrix(path0,"dump_Gc_","",str(i),False,False,1))
        Gf.append(load_matrix(path0,"dump_Gf_","",str(i),False,False,1))

        indBc = np.abs(Bc[i]).sum(axis=1)>0
        Bc_nonzRow.append( Bc[i][indBc,:])
#        Fc.append( np.dot(Bc_nonzRow[i], np.linalg.solve(K_reg[i],Bc_nonzRow[i].T)))
#        Lumped.append( np.dot(Bc_nonzRow[i], np.dot(K[i],Bc_nonzRow[i].T)))

        rhs.append(load_matrix(path0,"dump_rhs_","",str(i),False,False,1))
#        xx.append(load_matrix(path0,"dump_xxTest_","",str(i),False,False,1))
#        Kplus_f_test.append(load_matrix(path0,"dump_Kplus_f_test_","",str(i),False,False,1))


#        KplusBcT_p = BcKplus_List[i]

#        BcK_dense.append(load_matrix(path0,"dump_BcK_dense_","",str(i),False,False,1))
#        BcK_dense.append(np.dot(K[i],Bc_nonzRow[i].T).T)
        
        Gc.append(np.dot(Bc[i], R[i]))

        KplusBcT.append(load_matrix(path0,"dump_KplusBcT_","",str(i),False,False,1))
        KplusBcT_p.append(np.linalg.solve(K_reg[i],Bc_nonzRow[i].T))
#        BcKplus_tmp.append(np.linalg.solve(K_reg[i],Bc[i].T).T)

#        x_out.append(load_matrix(path0,"dump_x_out_","",str(i),False,False,1))
        
        Fc_p.append(np.dot(Bc_nonzRow[i],KplusBcT_p[i]))

#        iK_K = np.linalg.solve(K_reg[i],K[i])
#        K_iK_K = np.dot(K[i],iK_K)
#        del_ = np.linalg.norm(K_iK_K - K[i]  ) / np.linalg.norm(K[i])
#        print(del_)
#

        tmp_g = np.dot(Bc[i],np.linalg.solve(K_reg[i], rhs[i]))
        tmp_e = -np.dot(R[i].T,rhs[i])
        
        if (i == 0): 
            g_p  = tmp_g
            e_p = tmp_e;
        else:
            g_p += tmp_g;
            e_p = np.concatenate((e_p,tmp_e))


        print(' ...%d '%(i))

#    gc_p = np.concatenate((g_p,e_p)) 
#    gc_p = np.concatenate((gc_p,np.zeros(6))) 

    Gc_clust = load_matrix(path0,"dump_Gc_clust_","",str(0),False,False,1)
    Ac_clust = load_matrix(path0,"dump_Ac_clust_","",str(0),False,True,1)
    Fc_clust = load_matrix(path0,"dump_Fc_clust_","",str(0),False,True,1) 
    ker_GcTGc = load_matrix(path0,"dump_kerGc_","",str(0),False,False,1) 
#    gc = load_matrix(path0,"dump_gc_","",str(0),False,False,1)
#    lam_alpha = load_matrix(path0,"dump_lam_alpha_","",str(0),False,False,1)


#    lam_alpha_p = np.linalg.solve(Ac_clust, gc)

#    nLam = Bc[0].shape[0]
#    lam_p = lam_alpha_p[0:nLam]
##    alpha_p = lam_alpha[nLam:]
#    for i in range(nSub): 
#        print (" ! %d " % (i))
#        x10 = np.linalg.solve(K_reg[i],rhs[i])
#        x11 = np.linalg.solve(K_reg[i],np.dot(Bc[i].T,lam_p))
#
#        print alpha_p[(6*i):(6*(i+1))]
#        x2 = np.dot(R[i],alpha_p[(6*i):(6*(i+1))]) 
#
#        x_out_p.append(x10 - x11 + x2)

#        print( "||x_out - x_out_p || = %e " % np.linalg.norm(x_out[i] - x_out_p[i]))




    Ac_clust_python = np.hstack((Fc_clust,Gc_clust))

    Z = np.zeros((Gc_clust.shape[1],Ac_clust_python.shape[1]))
    print ( Z.shape)
    Ac_clust_python = np.vstack((Ac_clust_python,Z))

    Gf_clust = load_matrix(path0,"dump_Gf_clust_","",str(0),False,False,1)
#    test = load_matrix(path0,"dump_testXYZ_","",str(0),False,False,1)


#    KpOnes= load_matrix(path0,"dump_KplusONES_","",str(0),False,False,1)


#K_regD = K_reg[0]
#frhs = rhs[0]
#xxD = xx[0]
#RD = R[0]
#for i in range(1,nSub):
#    K_regD = block_diag(K_regD,K_reg[i]);
#    RD = block_diag(RD,R[i]);
#    frhs = np.concatenate((frhs,rhs[i]))
#    xxD = np.concatenate((xxD,xx[i]))
#

for i in range(nSub - 1):
    if (i == 0):
        Bc_g = np.hstack((Bc[0],Bc[1]))
    else:
        Bc_g = np.hstack((Bc_g,Bc[i+1]))

for i in range(nSub - 1):
    if (i == 0):
        Bf_g = np.hstack((Bf[0],Bf[1]))
    else:
        Bf_g = np.hstack((Bf_g,Bf[i+1]))

for i in range(nSub - 1):
    if (i == 0):
        Gf_g = Gf_p[0]+ Gf_p[1]
    else:
        Gf_g += Gf_p[i+1]




#Fc__ = np.dot(Bc_g,np.linalg.solve(K_regD,Bc_g.T))
#
#
#gc__ = np.dot(Bc_g,np.linalg.solve(K_regD,frhs))
#ec__ = - np.dot(RD.T,frhs)
#
#gc__ = np.concatenate((gc__,ec__))




#H = ker_GcTGc
#AA0 = np.hstack((Fc__,Gc_clust))
#AB1 = 
#
#
#ZZ1 = np.zeros((Gc_clust.shape[0], H.shape[1]))
#AA1 = np.vstack((ZZ1,H))
#AA01 = np.hstack((AA0,AA1))












#A0 = np.hstack((K_regD,Bc_g.T))
#
#nB = Bc_g.shape[0]
#Bc_Z = np.hstack((Bc_g,np.zeros((nB,nB))))
#
#crhs = np.zeros(nB);
#
#A = np.vstack((A0,Bc_Z))
#
#b = np.concatenate((frhs,crhs))
#
#x = np.linalg.solve(A,b)
#
#xxD = np.concatenate((xxD,crhs))




    
#Bc_g = np.hstack((Bc_g,Bc[2]))
#Bc_g = np.hstack((Bc_g,Bc[2])) 
#BcT_dense = load_matrix(path0,"dump_BcT_dense_","",str(0),True,True,1) 
#Fc_clust = load_matrix(path0,"dump_Fc_clust_","",str(0),True,True,1)
#Ac_clust = load_matrix(path0,"dump_Ac_clust_","",str(0),True,True,1)
#GcTGc = load_matrix(path0,"dump_GcTGc_clust_","",str(0),False,True,1) 
#GfTGf = load_matrix(path0,"dump_GfTGf_","",str(0),False,False,1) 
#iGfTGf = load_matrix(path0,"dump_iGfTGf_","",str(0),False,False,1) 
#ker_Ac = load_matrix(path0,"dump_ker_Ac_","",str(0),False,False,1) 
##KpBcT0 = load_matrix(path0,"dump_KplusBcT_","",str(0),False,False,1) 
##KpBcT1 = load_matrix(path0,"dump_KplusBcT_","",str(1),False,False,1) 
#
#
#dFc_eig = load_matrix(path0,"dump_Fc_clust_","",str(444),False,False,1)
##dFc_svd = load_matrix(path0,"dump_Fc_clust_","",str(555),False,False,1)
#dAc_eig = load_matrix(path0,"dump_Ac_clust_","",str(444),False,False,1)
##dAc_svd = load_matrix(path0,"dump_Ac_clust_","",str(555),False,False,1)
#
#
#GfTGf_  = np.zeros((GfTGf.shape[0],GfTGf.shape[0]))
#
#
#
#
#
#
#for d in range(nSub):
#    GfTGf_ += np.dot(Gf[d].T,Gf[d])
#
#
#
#
#if False:
#    plt.subplot(1,3,1)
#    if GcTGc.shape[0] < 100:
#        markersize_ = 3
#    else:
#        markersize_ = 0.7 
#    plt.spy(GcTGc, markersize=markersize_)
#    plt.xlabel("nnz = %d" % (GcTGc.nonzero()[0].shape[0]))
#    plt.subplot(1,3,2)
#    if Fc_clust.shape[0] < 100:
#        markersize = 3
#    else:
#        markersize = 0.7 
#    plt.spy(Fc_clust, markersize=markersize_)
#    plt.xlabel("nnz = %d" % (Fc_clust.nonzero()[0].shape[0]))
#    plt.subplot(1,3,3)
#    if Ac_clust.shape[0] < 100:
#        markersize_ = 3
#    else:          
#        markersize_ = 0.7 
#    plt.spy(Ac_clust, markersize=markersize_) 
#    plt.xlabel("nnz = %d" % (Ac_clust.nonzero()[0].shape[0]))
#    plt.show()
#
##Bc_from_Rt = []
##for i in range(1,14):
##    Bc_from_Rt.append( load_matrix(path0,"dump_Bc_from_Rt_","",str(i),False,False,1) )
##
#
## Gc_ = load_matrix(path0,"dump_Gc_i_","",str(0),False,False,1)
#
#
#
#
##BcT_dense = load_matrix(path0,"dump_BcT_dense_","",str(0),True,True,1) 
#
#
#K_test= []
#Kplus_K_test = []
#K_Kplus_K_test = []
#K_reg_test = []
#K_reg_SF = []
#x_test = []
#
#
#for i in range(4): 
#
#    K_test.append(load_matrix(path0,"dump_K_dense_","",str(i),False,True,1))
#    K_reg_test.append(load_matrix(path0,"dump_K_reg_","",str(i),False,True,1)) 
#    K_reg_SF.append(load_matrix(path0,"dump_K_reg_SF_","",str(i),False,True,1)) 
#    Kplus_K_test.append(load_matrix(path0,"dump_Kplus_K_","",str(i),False,False,1))
#    K_Kplus_K_test.append(load_matrix(path0,"dump_K_Kplus_K_","",str(i),False,False,1))
#
#    #KKpK = np.dot(K_test[i], np.linalg.solve(K_reg_test[i],K_test[i]))
#    KKpK = np.dot(K[i], np.linalg.solve(K_reg[i],K[i]))
#    print "norm = %3.8e \n" % np.linalg.norm(KKpK - K[i])
#


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##plt.spy(Fc_clust,markersize = .8);plt.show()
#
##Gc_ = load_matrix(path0,"dump_Gc_i_","",str(0),True,True,1)
#
#
#
##r = sparse.csgraph.reverse_cuthill_mckee(Ac_clust.tocsr(), symmetric_mode=True)
##Ac_clust = Ac_clust.toarray() 
###
##P,L,U= scipy.linalg.lu(Ac_clust)
##nnz0 = L.nonzero()[0].shape[0] +  U.nonzero()[0].shape[0]
##
##
###
###
##AcR = Ac_clust[np.ix_(r,r)]
##PR,LR,UR = scipy.linalg.lu(AcR)
##nnzR = LR.nonzero()[0].shape[0] +  UR.nonzero()[0].shape[0]
###
###
##plt.subplot(2,2,1)
##plt.spy(L,markersize=0.1);
##plt.subplot(2,2,2)
##plt.spy(U,markersize=0.1);
##plt.subplot(2,2,3)
##plt.spy(LR,markersize=0.1);
##plt.subplot(2,2,4)
##plt.spy(UR,markersize=0.1);
#
##print ("nnz = %d, nnz(reordered) = %d ") % (nnz0, nnzR) 
#
#
##plt.show()
#
##ker_Ac = load_matrix(path0,"dump_ker_Ac_","",str(0),False,True,1) 
##ker_GcTGc = load_matrix(path0,"dump_ker_GcTGc_","",str(0),False,True,1) 
##R0 = load_matrix(path0,"dump_R_","",str(0),False,True,1) 
#
##Gc_H = np.dot(GcTGc.toarray(),ker_GcTGc)
#
##r = sparse.csgraph.reverse_cuthill_mckee(Ac_clust.tocsr(), symmetric_mode=True)
##Ac = Ac_clust.toarray()[np.ix_(r,r)]
##plt.subplot(1,2,1)
##plt.spy(Ac_clust ,markersize = 2.0)
##plt.subplot(1,2,2)
##plt.spy(Ac,markersize = 0.125)
#
#
#
##Fc_python_List = []
#
##if 0:
##    Fc_clust_python = np.zeros((Bct_list[i].shape[0], Bct_list[i].shape[0]))
##    for i in range(nSub):
##        Bc = Bct_list[i].toarray()
##        indBc = np.abs(Bc).sum(axis=1)>0
##        Bc_red = Bc[indBc,:]
##        BcKplus = BcKplus_List[i]
##
##        Bf = Bf_List[i].toarray()
##        indBf = np.abs(Bf).sum(axis=1)>0
##        Bf_red = Bf[indBf,:]
##
##        Rc = RList[i].toarray()
##        
##        
##
##        if (i == 0):
##            Gf_clust_python = np.dot(Bf,Rc)
##            Gc_clust_python = np.dot(Bc,Rc)
##        else:
##            Gf_clust_python = np.hstack((Gf_clust_python,np.dot(Bf,Rc)))
##            Gc_clust_python = np.hstack((Gc_clust_python,np.dot(Bc,Rc)))
##        indBcKplus = np.abs(BcKplus).sum(axis=1)>0
##        BcKplus = BcKplus[indBcKplus,:] 
##        BcKplus_python = np.linalg.solve(K_reg_List[i],Bc_red.T)
##        BcKplus_ = np.linalg.solve(K_reg_List[i],Bc.T)
##        Fc_i = np.dot(Bc_red,BcKplus_python)
##        Fc_clust_python += np.dot(Bc,BcKplus_)
##        Fc_python_List.append(Fc_i)
##
##    for ii in range(nSub):
##     
##        ttt = Gc_List[ii][np.abs(Gc_List[ii]).sum(axis=1)>0,:] - GcList[ii]
##        print np.linalg.norm(ttt)
##
##
##    for ii in range(nSub):
##        ddd0 = np.linalg.norm(Fc_python_List[ii] - Fc_List[ii])
##        ddd1 = np.linalg.norm(Fc_python_List[ii])
##        print "|Fc_python - Fc_myAp|/|Fc_python|",ddd0 / ddd1  
##
##
##    Fc_clust =  load_matrix(path0,"dump_Fc_clust_","",0,False,True,1)
##    Gc_clust =  load_matrix(path0,"dump_Gc_clust_","",0,False,False,1)
##    Gf_clust =  load_matrix(path0,"dump_Gf_clust_","",0,False,False,1)
##    Ac_clust =  load_matrix(path0,"dump_Ac_clust_","",0,False,True,1)
##    Ac_clust_python = np.hstack((Fc_clust_python,Gc_clust_python))
##
##    Z = np.zeros((Gc_clust_python.shape[1],Ac_clust.shape[1]))
##    print ( Z.shape)
##    Ac_clust_python = np.vstack((Ac_clust_python,Z))
##
##
##    ddd0 = np.linalg.norm(Fc_clust - Fc_clust_python)
##    ddd1 = np.linalg.norm(Fc_clust)
##    print "|Fc_clust_python - Fc_clust_myAp|/|Fc_clust_python|",ddd0 / ddd1  
##
##    ddd0 = np.linalg.norm(Gc_clust - Gc_clust_python)
##    ddd1 = np.linalg.norm(Gc_clust)
##    print "|Gc_clust_python - Gc_clust_myAp|/|Gc_clust_python|",ddd0 / ddd1  
##    
##    ddd0 = np.linalg.norm(Gf_clust - Gf_clust_python)
##    ddd1 = np.linalg.norm(Gf_clust)
##    print "|Gf_clust_python - Gf_clust_myAp|/|Gf_clust_python|",ddd0 / ddd1  
##
##    ddd0 = np.linalg.norm(Ac_clust - Ac_clust_python)
##    ddd1 = np.linalg.norm(Ac_clust)
##    print "|Ac_clust_python - Ac_clust_myAp|/|Ac_clust_python|",ddd0 / ddd1  
##
##K = []
#
#
#
##plt.subplot(1,2,1)
##plt.spy(Gf_clust_python,markersize=1)
##plt.subplot(1,2,2)
##plt.spy(Gf_clust,markersize=1)
##plt.show()

