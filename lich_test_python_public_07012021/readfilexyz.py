import numpy as np
#import neighlisttmp as nblist
def readfile(filename):
    with open(filename, 'r') as f:
        # function [watoms, coord, latoms,coordlayer, xdim,ydim,zdim]=   readfile(filename)
        line1in = f.readline()
        line1 = line1in.split()
        print('Line 1 contents: ',line1)
        atoms = int(line1[0])
        print('Number of atoms',atoms)
        line2in = f.readline()
        line2=line2in.split()
        xdim  = float(line2[6]) # 7th column
        ydim  = float(line2[7]) # 8th column
        zdim  = float(line2[8]) # 9th column
        print('Cell dimensions:', xdim, ydim, zdim)
        coord=[]
        coordlayer=[]
        watoms=0
        latoms=0
        line  = f.readline()
        while len(line)>0:
            A = line.split()
            t=A[0]
            x=float(A[1])
            y=float(A[2])
            z=float(A[3])
            # Apply minimum image convention
            x = x - xdim * np.round(x / xdim)
            y = y - ydim * np.round(y / ydim)
            z = z - zdim * np.round(z / zdim)
            if (t=="mW") or (t=="OW"):
                coord.append([x, y, z])
                watoms=watoms+1
            else:
                coordlayer.append([x, y, z])
                latoms=latoms+1
            line  = f.readline()
        print('Water atoms:', watoms, 'Surface atoms:', latoms)
        print('Total atoms:', (watoms+latoms), 'N:', atoms)
        if atoms != (watoms+latoms):
            exit('Wrong number of atoms read in, exiting.')
        #print(coord)
        #print(coordlayer)
        return coord, coordlayer, xdim, ydim, zdim, watoms, latoms, line2in


#[coord, coordlayer, xdim, ydim, zdim, line2in] = readfile('coords.dat')
#print('From function readfile:', coord)
#[indexi, neigh] = nblist.neighbours(coord, xdim, ydim, zdim)
#print('Neighbour indices from nblist:\n', indexi)