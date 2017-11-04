from __future__ import print_function 
import numpy as np
from scipy import sparse
import myModul as mM
import config_espreso_python

#load_matrix0(path,str0,i,j,makeSparse,makeSymmetric,offset):


path0 = '/data_space/WorkSpace/htfeti_app/gitRepo/hfls/data/'

j="0"

#H = mM.load_matrix0(path0,"dump_ker_Ac_","",str(j),False,False,1)
H = mM.load_matrix0(path0,"dump_kerGc_","",str(j),False,False,1)
Ac = mM.load_matrix0(path0,"dump_Ac_clust_","",str(j),False,True,1)




HHt = np.dot(H,H.T)

Ac_rho = Ac.copy()


nA = Ac.shape[0]
nH = HHt.shape[0]
ind1 = np.arange(nA - nH, nA)

Ac_rho[np.ix_(ind1,ind1)] += HHt


#Sc = mM.load_matrix0(path0,"dump_S_new_","",str(j),False,False,1)
#Sc = Sc + Sc.T - np.diag(Sc.diagonal())
#Rc = mM.load_matrix0(path0,"dump_R_s_new_","",str(j),False,False,1)
#
#
#K = mM.load_matrix0(path0,"dump_K_new_","",str(j),False,True,1)
#R = mM.load_matrix0(path0,"dump_R_new_","",str(j),False,False,1)
#Rs = mM.load_matrix0(path0,"dump_R_s_new_","",str(j),False,False,1)
#Rr = mM.load_matrix0(path0,"dump_R_r_new_","",str(j),False,False,1)
#Krr = mM.load_matrix0(path0,"dump_K_rr_new_","",str(j),False,True,1)
#Kss = mM.load_matrix0(path0,"dump_K_ss_new_","",str(j),False,True,1)
#Krs = mM.load_matrix0(path0,"dump_K_rs_new_","",str(j),False,False,1)
#KrsRs = mM.load_matrix0(path0,"dump_K_rsRs_new_","",str(j),False,False,1)
##                                dump_K_rsRs_new_0.txt
#
#
#A11 = np.hstack([Krr,Krs])
#A22 = np.hstack([Krs.T,Kss])
#K_modif = np.vstack([A11,A22])
#
#
#Rr_ = np.linalg.solve(Krr,np.dot(Krs,Rs))


#R_modif = np.vstack([Rr,Rs])
