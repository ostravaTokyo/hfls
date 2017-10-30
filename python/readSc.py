
import numpy as np
from scipy import sparse
import myModul as mM
import config_espreso_python




path0 = '/data_space/WorkSpace/htfeti_app/gitRepo/hfls/data/'

j="0"
Sc = mM.load_matrix0(path0,"dump_Sc_clust_","",str(j),False,False,1)
