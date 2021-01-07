import numpy as np


# This version has a real cell list to speed up neighbour list generation

def neighbours(coord, xdim, ydim, zdim):
    N = len(coord)
    neigh_list = -1 * np.ones((N, 4), dtype='i')  # -1 marks empty neighbour slot
    atom_cell = np.zeros((N), dtype='i')  # Cell location for each atom
    print('Neighbourlist generation:')  # \n 0 %')
    rcut = 3.5
    ncellsx = int(np.floor(xdim / rcut))
    ncellsy = int(np.floor(ydim / rcut))
    ncellsz = int(np.floor(zdim / rcut))
    ncells = ncellsx * ncellsy * ncellsz
    print(f'Generating a cell list with {ncells} cells, x:{ncellsx}, y:{ncellsy}, z:{ncellsz}')

    # Find the correct cell for each atom
    print('Finding atomic positions in the cells')
    ilist = np.zeros((ncells), dtype='i')  # Current number of atoms in each cell
    for ni in range(0, N):
        # Find correct cell, to later add to that cell list cell
        #atomx = int(np.floor((coord[ni][0]) / xdim * ncellsx))  # For a cell from 0 to xdim
        #atomy = int(np.floor((coord[ni][1]) / ydim * ncellsy))
        #atomz = int(np.floor((coord[ni][2]) / zdim * ncellsz))
        atomx = int(np.floor((coord[ni][0]+0.5*xdim) / xdim * ncellsx))# +0.5*xdim for symmetric -xdim/2 to xdim2/2 cell
        atomy = int(np.floor((coord[ni][1] + 0.5 * ydim) / ydim * ncellsy))
        atomz = int(np.floor((coord[ni][2] + 0.5 * zdim) / zdim * ncellsz))
        atom_cell[ni] = atomx*ncellsy*ncellsz + atomy*ncellsz + atomz # Cell number
        ilist[atom_cell[ni]] += 1
        #print('Atom_cell for atom', ni, 'is',atom_cell[ni])
        #print('Coords:', coord[ni][0], coord[ni][1], coord[ni][2])
        #print('atomx, y, z:', atomx, atomy, atomz)
    print('Atom cells:', atom_cell)

    maxatcell = np.max(ilist)
    print('Maximum atoms in a cell:', maxatcell)
    cell_list = -1 * np.ones((ncells, maxatcell+1), dtype='i')  # -1 marks empty slot. This has local atoms
    cell_list_neigh = -1 * np.ones((ncells, 27*maxatcell), dtype='i')  # This has also neighbour cell atoms
    cell_location = np.zeros((ncells, 3), dtype='i')  # Keeps x, y, z order of cells

    #   - Preparation: -
    # Generate the x, y, z list of cells (for example, in a 15**3 list you could have (7, 7, 7) around origo.)
    # Using cell location, update cell_neigh to have the surrounding 26 cells for each cell, using cell_location
    #  not forgetting the periodic boundaries for the first and last location in each dimension.
    #  Simple if nx>ncells, nx = 0 etc.
    # Give all atoms their cell number based on coord, into atom_cell
    ci = 0  # index of the current cell, for setting the location counters
    for cx in range(0, ncellsx):
        for cy in range(0, ncellsy):
            for cz in range(0, ncellsz):
                cell_location[ci, 0] = cx
                cell_location[ci, 1] = cy
                cell_location[ci, 2] = cz
                ci += 1

    print('Setting up cell neighbours')
    cell_neigh = np.zeros((ncells, 26), dtype='i') # Set the cell_neigh array with all the neighbour cell indices
    for ci in range(0, ncells):
        cx = cell_location[ci, 0]  # These have current cell's position
        cy = cell_location[ci, 1]
        cz = cell_location[ci, 2]
        czmodminus = ((cz - 1) % ncellsz) - cz  # There pre-calculate the modulo parts needed in all the following terms
        czmodplus  = ((cz + 1) % ncellsz) - cz
        cymodminus = (((cy - 1) % ncellsy) - cy) * ncellsz
        cymodplus  = (((cy + 1) % ncellsy) - cy) * ncellsz
        cxmodminus = (((cx - 1) % ncellsx) - cx) * ncellsy * ncellsz
        cxmodplus  = (((cx + 1) % ncellsx) - cx) * ncellsy * ncellsz
        
        cell_neigh[ci, 0] = ci + czmodminus  # -z cell. Modulo needed for periodic boundaries
        cell_neigh[ci, 1] = ci + czmodplus  # +z cell.
        cell_neigh[ci, 2] = ci + cymodminus  # -y cell.
        cell_neigh[ci, 3] = ci + cymodplus  # +y cell.
        cell_neigh[ci, 4] = ci + cxmodminus  # -x cell.
        cell_neigh[ci, 5] = ci + cxmodplus  # +x cell.
        cell_neigh[ci, 6] = ci + cymodminus + czmodminus  # -y, -z cell.
        cell_neigh[ci, 7] = ci + cymodminus + czmodplus  # -y, +z cell.
        cell_neigh[ci, 8] = ci + cymodplus + czmodminus  # +y, -z cell.
        cell_neigh[ci, 9] = ci + cymodplus + czmodplus  # +y, +z cell.
        cell_neigh[ci, 10] = ci + cxmodminus + czmodminus  # -x, -z cell.
        cell_neigh[ci, 11] = ci + cxmodminus + czmodplus  # -x, +z cell.
        cell_neigh[ci, 12] = ci + cxmodplus + czmodminus  # +x, -z cell.
        cell_neigh[ci, 13] = ci + cxmodplus + czmodplus  # +x, +z cell.
        cell_neigh[ci, 14] = ci + cxmodminus + cymodminus  # -x, -y cell.
        cell_neigh[ci, 15] = ci + cxmodminus + cymodplus   # -x, +y cell.
        cell_neigh[ci, 16] = ci + cxmodplus + cymodminus  # +x, -y cell.
        cell_neigh[ci, 17] = ci + cxmodplus + cymodplus   # +x, +y cell.
        cell_neigh[ci, 18] = ci + cxmodminus + cymodminus + czmodminus  # -x, -y, -z cell.
        cell_neigh[ci, 19] = ci + cxmodminus + cymodplus + czmodminus  # -x, +y, -z cell.
        cell_neigh[ci, 20] = ci + cxmodplus + cymodminus + czmodminus  # +x, -y, -z cell.
        cell_neigh[ci, 21] = ci + cxmodplus + cymodplus + czmodminus  # +x, +y, -z cell.
        cell_neigh[ci, 22] = ci + cxmodminus + cymodminus + czmodplus  # -x, -y, +z cell.
        cell_neigh[ci, 23] = ci + cxmodminus + cymodplus + czmodplus  # -x, +y, +z cell.
        cell_neigh[ci, 24] = ci + cxmodplus + cymodminus + czmodplus  # +x, -y, +z cell.
        cell_neigh[ci, 25] = ci + cxmodplus + cymodplus + czmodplus  # +x, +y, +z cell.

    # Then setup local cell list
    print('Writing atom indices to the cell list')
    nlist = np.zeros((ncells), dtype='i')  # Current position to fill in cell_list column
    for ni in range(0, N):
        # Add atoms to the correct cell of the cell list
        cell_list[atom_cell[ni], nlist[atom_cell[ni]]] = ni   # Add the atom index to the correct cell list position
        nlist[atom_cell[ni]] += 1    # Update the next position to write to the cell list
    # Next another for-loop over cells cin, which pics the indices of neighbouring cells, and copies their atoms
    # to the neighbour lists
    #print('Writing atom indices to the extended cell list with neighbours')
    #nlistn = np.zeros((ncells), dtype='i')  # Current position to fill in cell_list_neigh column
    #for ci in range(0, ncells):  # Set the cell_neigh array with all the
    #    ineigh = 0
    #    while cell_list[ci, ineigh] > -1:  # Here copy the local atoms to cell_list_neigh
    #        cell_list_neigh[ci, ineigh] = cell_list[ci, ineigh]
    #        ineigh += 1
    #        nlistn[ci] += 1    # Update the number of items in the extended cell list
    #    for cin in range(0, 26):  # Go through all 26 neighbouring cells
    #        ineigh2 = 0  # Atom indices within each neighbouring cell list
    #        while cell_list[cell_neigh[ci, cin], ineigh2] > -1:  # Maybe like this...
    #            # Add atoms to cell_list_neigh
    #            # Update ineigh etc.
    #            cell_list_neigh[ci, ineigh] = cell_list[cell_neigh[ci, cin], ineigh2]
    #            ineigh2 += 1
    #            ineigh += 1
    #            nlistn[ci] += 1    # Update the number of items in the extended cell list

        #   - Generation of the final cell list: cell_list_neigh -
        # Add all atoms into their own cell list: cell_list
        # Then with a new for loop, check which cell the atom i is in from atom_cell
        # Then check from cell_neigh which are the neighbour cells
        # Actually just one line copying the atom list from these two to cell_list_neigh is enough
        #   - Generation of the neighbour list -
        # For each atom, copy atom from own column of cell_list_neigh,
        # then copy from indices found from cell neigh all atoms.

    #print(cell_list)
    #print(cell_list[2083])
    # indexi = []
    #neigh = [] NOT IN USE
    iprog = 0  # Progress counter
    # Find neighbours from the cell list
    print('Generating the neighbour list from the cell list')
    neigh_num = np.zeros((N), dtype='i')  # Number of neighbours for each atom
    neightot = 0   # Total number of neighbours found  
    for i in range(0, N):
        if (i + 1) % int(float(N) / 10) == 1:  # If 10 %, 20 %, ...
            print('Progress:', iprog, '%')
            iprog = iprog + 10
        ri = []  # List of neigbouring radii
        nneigh = 0  # Gathers the number of possible neighbours found
        maxr = 0  # Gathers the maximum distance in all cell_lists
        xi = coord[i][0]
        yi = coord[i][1]
        zi = coord[i][2]
        # First the local cell:
        for jc in range(0, nlist[atom_cell[i]]):  # Go through all local neighbours of atom i
            j = cell_list[atom_cell[i], jc]
            if cell_list[atom_cell[i], jc] == -1:  # End of filled list, not needed if nlist correct
                print('breaking local for loop for jc:', jc)
                break
            if (i != j):               
                dx = coord[j][0] - xi  # Missing the if statements
                dy = coord[j][1] - yi  # below which speed up
                dz = coord[j][2] - zi
                # Apply minimum image convention
                dxm = dx - xdim * np.round(dx / xdim)
                dym = dy - ydim * np.round(dy / ydim)
                dzm = dz - zdim * np.round(dz / zdim)
                # Norm of minimum image vector
                rsq = dxm * dxm + dym * dym + dzm * dzm
                r = rsq ** 0.5
                if (r > maxr):
                    maxr = r
                nneigh += 1
                if abs(r) < rcut:  # Orig. 3.6, why not rcut?
                    ri.append([r, j])
        # Next the neighbouring cells:
        ci = atom_cell[i] # Current cell
        for cin in range(0, 26):  # Go through all 26 neighbouring cells
            ineigh2 = 0  # Atom indices within each neighbouring cell list   #while cell_list[cell_neigh[ci, cin], ineigh2] > -1:  # Maybe like this...
            for jc in range(0, nlist[cell_neigh[ci, cin]]):  # Go through all members if this cell
                j = cell_list[cell_neigh[ci, cin], jc]
                if cell_list[cell_neigh[ci, cin], jc] == -1:  # End of filled list, not needed if nlist correct
                    print('breaking surrounding cells for loop for jc:', jc)
                    break                          
                dx = coord[j][0] - xi  # Missing the if statements
                dy = coord[j][1] - yi  # below which speed up
                dz = coord[j][2] - zi
                # Apply minimum image convention
                dxm = dx - xdim * np.round(dx / xdim)
                dym = dy - ydim * np.round(dy / ydim)
                dzm = dz - zdim * np.round(dz / zdim)
                # Norm of minimum image vector
                rsq = dxm * dxm + dym * dym + dzm * dzm
                r = rsq ** 0.5
                if (r > maxr):
                    maxr = r
                nneigh += 1
                if abs(r) < rcut:  # Orig. 3.6, why not rcut?
                    ri.append([r, j])

                    
        ri.sort()
        ri = ri[:4]
        # print(ri)
        # ind = []
        j = 0
        for ineigh in ri:  # This construct to make sure only running to the max number of neighbours
            neigh_list[i, j] = ineigh[1]
            # neigh.append(coord(ineigh[1]))
            j += 1
            neightot += 1
            neigh_num[i] += 1 # Update the number of neighbours for atom i
        #print('Max r of cell_list_neigh:', maxr, 'Number of possible neighbours:', nneigh )  # Debugging
        # print(ind)
        # indexi.append(ind)
        # for item in ind:
        #    neigh.append(coord[item])
        # neigh_list = indexi

    # print('Indices for neighbours:', neigh_list, 'Neighbour coords:')
    # print(neigh)
    avgneighn = neightot / N   # Average number of neighbours per atom
    print('Avg. # of neighbours:', avgneighn)
    print('Type and size of neigh_list:', type(neigh_list), neigh_list.shape)
    #print(neigh_list[:10])
    # print('Progress: 100 %')
    return neigh_list, N, neigh_num  # Note that names should now match the call!
