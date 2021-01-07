import numpy as np
#import neighlisttmp as nblist
def readfile(filename):
    typeow = 'OW'  # Marker for the correct atom type to read in
    with open(filename, 'r') as f:  # First read in just the cell size
        lines = f.read().splitlines()
        last_line = lines[-1]
    print(last_line)
    linecell = last_line.split()
    xdim  = 10*float(linecell[0])  # .gro has nm as units
    ydim  = 10*float(linecell[0]) 
    zdim  = 10*float(linecell[0])
    coord = [] 
    coordlayer = []
    watoms = 0
    latoms = 0
    with open(filename, 'r') as f:
        # function [watoms, coord, latoms,coordlayer, xdim,ydim,zdim]=   readfile(filename)
        line1in = f.readline()
        line1 = line1in.split()
        print('Line 1 contents: ',line1)
        line2in = f.readline()
        line2=line2in.split()
        print('Line 2 contents: ',line2)        
        atoms = int(line2[0])
        print('Number of atoms',atoms)        
        line  = f.readline()
        for i in range(atoms):
            if typeow in line:
                A = line.split()
                t=A[1]
                x=10*float(A[2])-0.5*xdim
                y=10*float(A[3])-0.5*ydim
                z=10*float(A[4])-0.5*zdim
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
        print('Water OW atoms:', watoms, 'Surface atoms:', latoms)
        print('Total atoms:', (watoms+latoms), 'Ntot (4 per water):', atoms)
        #if atoms != (watoms+latoms):
        #    exit('Wrong number of atoms read in, exiting.')

       # # Apply minimum image convention
       # for i in range(atoms)
       #     x = x - xdim * np.round(x / xdim)
       #     y = y - ydim * np.round(y / ydim)
       #     z = z - zdim * np.round(z / zdim)
        #print(coord)
        #print(coordlayer)
        return coord, coordlayer, xdim, ydim, zdim, watoms, latoms, line2in


#[coord, coordlayer, xdim, ydim, zdim, watoms, latoms, line2in] = readfile('initial.gro')
#print('From function readfile:', coord)
#[indexi, neigh] = nblist.neighbours(coord, xdim, ydim, zdim)
#print('Neighbour indices from nblist:\n', indexi)
