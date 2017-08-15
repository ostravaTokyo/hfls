import numpy as np
from scipy.linalg import block_diag


Kloc = np.array([[1.,-1.],[-1.,1.]])

Kploc = np.zeros((2,2))

Kploc[0,0] = 1;

K = block_diag(Kloc,Kloc)
Kp = block_diag(Kploc,Kploc)




KKpK = np.dot(K, np.dot(Kp,K))

print KKpK

B = np.array([[0,1.0,-1.0,0.]])

R = np.zeros((4,2))
R[:2,0] = 1.
R[2:,1] = 1.

G = np.dot(B,R)

F = np.dot(B, np.dot(Kp,B.T))

A_00 = np.hstack((F,G))

Z = np.zeros((1))

A_10 = np.hstack((G.T,Z))


A = np.vstack((A_00,A_10))
