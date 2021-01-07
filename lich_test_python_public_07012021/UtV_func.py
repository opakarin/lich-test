#Matlab function [Tx,disconnection_flag] = UtV_func(U,V)
import numpy as np
def UtV_func(U, V, i):

    #A = D[i,:,:]
    #B = D[j,:,:]
    U = U / ((sum(U ** 2) + np.finfo(float).eps) ** 0.5)
    V = V / ((sum(V ** 2) + np.finfo(float).eps) ** 0.5)

    Tx = np.matmul(np.transpose(U), V) # Previously X <- A, B
    #print('Early version of Tx: \n', Tx)
    indsij = np.where(abs(Tx + 1) < 1e-7)
    #if i < 15:
    #    print('Early version of Tx: \n', Tx)
    #    print('indsij of i:', i, indsij[:])
    ii = indsij[0]
    jj = indsij[1]
    #print('ii, jj, len: \n',ii, jj, len(ii), len(jj))
    #print('argwhere_ii, argwhere_jj:',np.argwhere(ii), np.argwhere(jj))
    if len(ii)<1:              # if there is no O-O bond
        disconnection_flag = 1 # there is no tetrahedral structure (disconnection in O-O)
    else:
        Tx = np.delete(Tx, ii[0], 0)  # X[ii(1),:]=[]
        Tx = np.delete(Tx, jj[0], 1)  # X[:,jj(1)]=[]
        disconnection_flag = 0
        #print('Tx after row and column deletion: \n', Tx)

    return Tx, disconnection_flag

