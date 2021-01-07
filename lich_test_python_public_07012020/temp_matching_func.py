# Matlab function [S, E] = temp_matching_func(neigh_dirs,neigh_list,neigh_num,atoms_num,score_min)
from score_func import score_func 
from UtV_func import UtV_func 
import numpy as np
import time

def temp_matching_func(neigh_dir,neigh_list,neigh_num,N,min_score):

 
    S = np.zeros((N, 4))
    E = np.zeros((N, 4))
    timefutv = 0.0  # Time counter for UtV_func
    timeconf = 0.0  # Time counter for score_func
    Ts = np.array([[-1, 0.5, 0.5],
            [0.5, -1, 0.5],
            [0.5, 0.5, -1]])
    Te = np.array([[0.5, -0.5, -0.5], 
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5]])
    P = np.array([[0, 0, 1, 0, 1, 0, 1, 0, 0],
        [  0, 0, 1, 1, 0, 0, 0, 1, 0],
        [  0, 1, 0, 0, 0, 1, 1, 0, 0],
        [  0, 1, 0, 1, 0, 0, 0, 0, 1],
        [  1, 0, 0, 0, 0, 1, 0, 1, 0],
        [  1, 0, 0, 0, 1, 0, 0, 0, 1]])
    for i in range(N):
        for j in range(4):
            jth_neigh_i = neigh_list[i, j]
            if jth_neigh_i > i:
                tic = time.perf_counter()               
                [Tx, discon_flag] = UtV_func(neigh_dir[i,:,:], neigh_dir[jth_neigh_i,:,:], i)
                toc = time.perf_counter()
                timefutv += (toc-tic)
                nU = neigh_num[i] - 1
                nV = neigh_num[jth_neigh_i] - 1
                # print('X, i, j: \n', X, i, j) # For debugging
                if discon_flag == 0:
                    tic = time.perf_counter()
                    [sscr, escr] = score_func(Tx, Ts, Te, nU, nV, P) 
                    toc = time.perf_counter()
                    timeconf += (toc-tic)
                    #print('Scoring atom: ', i)
                else:
                    neigh_list[i, j] = -1  # Not anymore counted as neighbours ## Should neigh_num be updated too? ##
                    sscr = 0
                    escr = 0
                S[i, j] = sscr
                E[i, j] = escr
                for k in range(4):
                    if neigh_list[jth_neigh_i, k] == i:
                        S[jth_neigh_i, k] = sscr
                        E[jth_neigh_i, k] = escr
#                        S[jth_neigh_i, neigh_list[jth_neigh_i, k]] = sscr
#                        E[jth_neigh_i, neigh_list[jth_neigh_i, k]] = escr

    E = E*((E >= S).astype(int)) # E[E<S]=0
    S = S*((S >= E).astype(int)) # S[S<E]=0
    E = E*((E >= min_score).astype(int)) # E[E<min_score]=0
    S = S*((S >= min_score).astype(int)) # S[S<min_score]=0
    
    # print('S[i, :] later: ',S[:10, :])
    # print('E[i, :] later: ',E[:10, :])            
    print ('Time spent on UtV_func, score_func:', timefutv, timeconf)            
    return S, E # Changed name from Sscore, Escore
