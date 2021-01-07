# Matlab function Labels = ice_labelling_func(Stg_score, Ecl_score, neigh_list, neigh_num, inds_ext, watoms_num)
# Version 3.1 with interfacial ice subtypes
import numpy as np
import time

import numpy as np


def ice_labelling_func(Stg_score, Ecl_score, neigh_list, neigh_num, N):
    
    labels = np.zeros((N), dtype='i')
    Stg_score = np.around(Stg_score+0.0001) # Sets to logical (0/1) values, +0.0001 needed for rounding 0.5 up.
    Ecl_score = np.around(Ecl_score+0.0001) # Sets to logical (0/1) values, +0.0001 needed for rounding 0.5 up.
    for i in range(N):
        if sum(Stg_score[i, :]) == 4 and sum(Ecl_score[i, :]) == 0: 
            labels[i] = 1  # Cubic ice
        elif sum(Stg_score[i, :]) == 3 and sum(Ecl_score[i, :]) == 1: 
            labels[i] = 2  # Hexagonal ice

    neighlab = np.zeros((N, 4), dtype='i')         
    for i in range(N):
        if sum(Stg_score[i, :]) + sum(Ecl_score[i, :]) == 0 or labels[i] > 0:
            continue    # Only deal with potential interfacial ice here
        for j in range(neigh_num[i]):
            neighlab[i, j] = labels[neigh_list[i, j]]
        if any(l == 1 for l in neighlab[i, :]) and any(l == 2 for l in neighlab[i, :]):
            labels[i] = 3  # Mixed
        elif any(l == 1 for l in neighlab[i, :]):  
            labels[i] = 4  # CI
        elif any(l == 2 for l in neighlab[i, :]):  
            labels[i] = 5  # HI

    for i in range(N):
            if sum(Ecl_score[i, :]) >= 3:  # CH
                   labels[i] = 6    

    neighlab = np.zeros((N, 4), dtype='i')
    for i in range(N):
        if labels[i] > 0:
            continue     # Only deal with potential interfacial ice here
        iint = False
        for j in range(neigh_num[i]):
            neighlab[i, j] = labels[neigh_list[i, j]]
            if neighlab[i, j] == 3 or neighlab[i, j] == 4 or neighlab[i, j] == 5 or neighlab[i, j] == 6:  # one of neighbours is interfacial or CH
                if Stg_score[i, j] + Ecl_score[i, j] > 0:  # and staggered or eclipsed configuration to that direction
                    iint = True
        if iint == True:
            labels[i] = 7  # I   
            
    return labels
