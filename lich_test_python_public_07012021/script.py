# Main script for running IcePro analysis in Python3
import numpy as np
import time
#from readfilexyz_outputbox import readfile
from readfilexyz import readfile
from neighlistcell import neighbours
from temp_matching_func import temp_matching_func
from find_neigh_dir import find_neigh_dir
from ice_labelling_func import ice_labelling_func
#filename = 'initial.gro.xyz'
filename = 'tmp_52000sym.xyz'
#filename = 'ice_c.xyz'
# Read in the coordinate file
tic0 = time.perf_counter() # Starting time
[coord, coordlayer, xdim, ydim, zdim, watoms, latoms, line2in] = readfile(filename)
# Generate neighbour list for water oxygen atoms
tic = time.perf_counter()
[neigh_list, N, neigh_num] = neighbours(coord, xdim, ydim, zdim)
toc = time.perf_counter()
timenlist = toc - tic # Time counter for neighbours()

tic = time.perf_counter()
neigh_dir = find_neigh_dir(coord, neigh_list, xdim, ydim, zdim, N)
toc = time.perf_counter()
timendir = toc - tic # Time counter for find_neigh_dir

print('Giving scores to configurations')
tic = time.perf_counter()
min_score = 0.5; # the least considerable similarity score between the structure of atoms in the box and a template
[Stg_score, Ecl_score] = temp_matching_func(neigh_dir,neigh_list,neigh_num,N,min_score)
toc = time.perf_counter()
timescores = toc - tic # Time counter for SCORES()

tic = time.perf_counter()
print('Calculating ice types from configuration scores')
labels = ice_labelling_func(Stg_score, Ecl_score, neigh_list, neigh_num, N)  # Labels ice types based on scores

# ncubic = 0  # These are for the outputbox version
# nhex = 0
# nint = 0
# nliq = 0
# for i in range(N):
#     if (coord[i][0] > outputbox[0] and coord[i][0] < outputbox[1] \
#         and coord[i][1] > outputbox[2] and coord[i][1] < outputbox[3] \
#         and coord[i][2] > outputbox[4] and coord[i][2] < outputbox[5]):
#         if labels[i] == 1:
#             ncubic += 1
#         elif labels[i] == 2:
#             nhex += 1
#         elif labels[i] == 3:
#             nint += 1
#         else:
#             nliq += 1

ncubic = len(np.argwhere(labels[:] == 1))  # These for the normal, full cell version
nhex = len(np.argwhere(labels[:] == 2))
nmix = len(np.argwhere(labels[:] == 3))
nci = len(np.argwhere(labels[:] == 4))
nhi = len(np.argwhere(labels[:] == 5))
nch = len(np.argwhere(labels[:] == 6))
nint = len(np.argwhere(labels[:] == 7))

nliq = N - ncubic - nhex - nmix - nci - nhi - nch - nint  # If adding other atom types, change accordingly   

toc = time.perf_counter()
timelabel = toc - tic # Time counter for ice type labeling
        
file_out = filename + '_ice'
print('Writing out the analyzed structure to', file_out)
with open(file_out, 'w') as fw:
    fw.write(str(N)+' \n')
    fw.write(line2in)
    for i in range(N):
        if labels[i]==1:
            fw.write('C '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        elif labels[i]==2:
            fw.write('H '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        elif labels[i]==3:
            fw.write('IM '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        elif labels[i]==4:
            fw.write('IC '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        elif labels[i]==5:
            fw.write('IH '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        elif labels[i]==6:
            fw.write('CH '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        elif labels[i]==7:
            fw.write('In '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
        else:
            fw.write('L '+ str(coord[i][0])+' '+ str(coord[i][1])+' '+ str(coord[i][2])+' \n')
    #for i in range(layeratoms):
    #   fw.write('X'+ str(coordlayer[i][0])+ str(coordlayer[i][1])+ str(coordlayer[i][2])+' \n')

toc = time.perf_counter()
timetot = toc - tic0 # Total time
file_dat = filename + '.dat'
print('Writing out the number of cubic, hex, mixed int., cubic int., hex. int., clathrate, interfacial, liquid to', file_dat)
print(ncubic, nhex, nmix, nci, nhi, nch, nint, nliq)
with open(file_dat, 'w') as fw:
    fw.write(str(ncubic) +' ' + str(nhex) +' ' + str(nmix) +' ' + str(nci) +' ' + str(nhi) +' ' + str(nch) +' ' + str(nint) +' ' + str(nliq) +' \n')
    
print('Time used for subroutines:')
print('Neighbour list:', timenlist)
print('Neighbour directions:', timendir)
print('Configuration scoring:', timescores)
print('Ice type labeling:', timelabel)
print('Total time:', timetot)
