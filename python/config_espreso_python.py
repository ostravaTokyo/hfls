#config

eps_dual_feti   = 1e-13
maxIt_dual_feti = 500


iterative_Kplus                     = False# True, False



eps_iter_Kplus                      = 1e-15
single_precision                    = False# True, False
methodToImproveSolByIterMethod      = 'pcg_x'   # 'cg_x', 'pcg_x',       
precondFrom_Areg_orA                = False
precondPrimalSystem                 = 'none'    # 'diag', 'LU_SP', 'ILU,  'none'

#       WARNING:  ILU (incomplete LU) can generate 
#       unsymmetric preconditioner matrix!!!


# if 'cg_x' or is set
mult_Areg_or_A_RRt                  = True








###############################################################################
###############################################################################
###############################################################################
###############################################################################
if single_precision==True:
    iterative_Kplus==True    
if precondPrimalSystem ==  'LU_SP':
    single_precision = True

#if methodToImproveSolByIterMethod = 'pcg_x':
#    precondFrom_Areg_orA = False
    
    
    
# methodToImproveSolByIterMethod      = 'cg_x'
# precondFrom_A_or_Areg               = True
# mult_Areg_or_A_RRt                  = False
#
#
#
#
    
#'cg_dx''cg_dx' 
