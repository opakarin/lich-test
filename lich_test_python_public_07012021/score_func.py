import numpy as np
from scipy.spatial import distance


def score_func(Tx, Ts, Te, nU, nV, P):

    lamb = 0.15
    DS = np.zeros((3, 3))
    DE = np.zeros((3, 3))

    if (nV > 0 and nU > 0):
        arr1 = np.transpose(Tx[:nU, :nV])
        arr2 = np.transpose(Ts[:nU, :])
        #print('Sizes of Tx_prime and Ts, nU, nV:', arr1.shape, arr2.shape, nU, nV)
        #print('Sizes of Tx and Ts_prime, nU, nV:', arr1.shape, arr2.shape, nU, nV)
        DStmp = distance.cdist(np.transpose(Tx[:nU, :nV]), np.transpose(Ts[:nU, :]))**2
        DEtmp = distance.cdist(np.transpose(Tx[:nU, :nV]), np.transpose(Te[:nU, :]))**2
        #DStmp = distance.cdist(Tx[:nU, :nV], np.transpose(Ts[:nV, :nU]))**2 # Only testing
        #print('DStmp:',  DStmp)
        DS[:nV,:] = DStmp
#        DS[:nV,:] = distance.cdist(np.transpose(Tx[:nU, :nV]), Ts[:nV, :])**2
        #DEtmp = distance.cdist(np.transpose(Tx[:nU, :nV]), np.transpose(Te[:nV, :nU]))**2 # Only testing
        #print('DStmp:',  DEtmp)
        DE[:nV,:] = DEtmp
#        DE[:nV,:] = distance.cdist(np.transpose(Tx[:nU, :nV]), Te[:nV, :])**2
        DS1D = np.ravel(DS, order='F')
        DS1D = DS1D.reshape(-1,1)
        DE1D = np.ravel(DE, order='F')
        DE1D = DE1D.reshape(-1,1)
        dsmin = np.amin(np.matmul(P, DS1D))/(lamb*nV*nU)
        demin = np.amin(np.matmul(P, DE1D))/(lamb*nV*nU)
        sscr = np.exp(-dsmin)
        escr = np.exp(-demin)
    else:
        sscr=0
        escr=0

 
    return sscr, escr
