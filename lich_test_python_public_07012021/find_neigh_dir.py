#function Neigh_dir = find_neigh_dir(coord,Neigh_list)
import numpy as np

def find_neigh_dir(coord, neigh_list, xdim, ydim, zdim, N):

    print('Calculating neighbour directions')
    neigh_dir = np.zeros((N, 3, 4))  # 3d tensor of vectors to 4 neighbours
    for i in range(N):
        xi = coord[i][0]
        yi = coord[i][1]
        zi = coord[i][2]
        #print('i:', i)
        #print(type(neigh_list))
#        print(neigh_list[1:10,:])
        #print('Neigh list columns', len(neigh_list[i]))
        #for j in range(len(neigh_list[i,:])):
        for j in range(4):
            if neigh_list[i][j] != -1:
                jind = neigh_list[i,j]
                dx = coord[jind][0] - xi
                dy = coord[jind][1] - yi
                dz = coord[jind][2] - zi
                # Apply minimum image convention
                dxm = dx - xdim * np.round(dx / xdim)
                dym = dy - ydim * np.round(dy / ydim)
                dzm = dz - zdim * np.round(dz / zdim)
                neigh_dir[i, 0 ,j] = dxm
                neigh_dir[i, 1, j] = dym
                neigh_dir[i, 2, j] = dzm
    return neigh_dir
