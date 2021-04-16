function [neigh_list, neigh_dir, neigh_num, inds_ext] = neigh_info_func(coordinates,box_x,box_y,box_z,cutoff)
% This function returns neighbour's information for water oxygen atoms
% [list_of_neighbours, direction_of_bonds, num_of_neighbours, 
% indexes_of_neighbouring_atoms_repeated_in_periodic_images]...
% = neigh_info_func(coordinates,box_size1,box_size2,box_size3,cutoff)

%% finding indexes of atoms which their periodic images are in the neighbourhood of the box (dis. less than cutoff)
inds_x_righ = find(coordinates(:,1)> box_x-cutoff);
inds_x_left = find(coordinates(:,1)< cutoff);
inds_y_righ = find(coordinates(:,2)> box_y-cutoff);
inds_y_left = find(coordinates(:,2)< cutoff);
inds_z_righ = find(coordinates(:,3)> box_z-cutoff);
inds_z_left = find(coordinates(:,3)< cutoff);

%%% corners
inds_xy_rightop = find(coordinates(:,1)> box_x-cutoff & coordinates(:,2)> box_y-cutoff);
inds_xy_righbot = find(coordinates(:,1)> box_x-cutoff & coordinates(:,2)< cutoff);
inds_xy_lefttop = find(coordinates(:,1)< cutoff & coordinates(:,2)> box_y-cutoff);
inds_xy_leftbot = find(coordinates(:,1)< cutoff & coordinates(:,2)< cutoff);

inds_xz_rightop = find(coordinates(:,1)> box_x-cutoff & coordinates(:,3)> box_z-cutoff);
inds_xz_righbot = find(coordinates(:,1)> box_x-cutoff & coordinates(:,3)< cutoff);
inds_xz_lefttop = find(coordinates(:,1)< cutoff & coordinates(:,3)> box_z-cutoff);
inds_xz_leftbot = find(coordinates(:,1)< cutoff & coordinates(:,3)< cutoff);

inds_yz_rightop = find(coordinates(:,2)> box_y-cutoff & coordinates(:,3)> box_z-cutoff);
inds_yz_righbot = find(coordinates(:,2)> box_y-cutoff & coordinates(:,3)< cutoff);
inds_yz_lefttop = find(coordinates(:,2)< cutoff & coordinates(:,3)> box_z-cutoff);
inds_yz_leftbot = find(coordinates(:,2)< cutoff & coordinates(:,3)< cutoff);
%%%

%inds_ext = [inds_x_righ;inds_x_left;inds_y_righ;inds_y_left;inds_z_righ;inds_z_left]; % the indexes of all repeating ...
% ... atoms in the extended box

%%% corners
inds_ext = [inds_x_righ;inds_x_left;inds_y_righ;inds_y_left;inds_z_righ;inds_z_left;...
    inds_xy_rightop; inds_xy_righbot; inds_xy_lefttop; inds_xy_leftbot; inds_xz_rightop; inds_xz_righbot; ...
    inds_xz_lefttop; inds_xz_leftbot; inds_yz_rightop; inds_yz_righbot; inds_yz_lefttop; inds_yz_leftbot ]; % the indexes of all repeating ...
%%%

%% extending the box in order to consider the periodic neighborhoods
X_shift_rl = [coordinates(inds_x_righ,1)-box_x coordinates(inds_x_righ,2:3)];
X_shift_lr = [coordinates(inds_x_left,1)+box_x coordinates(inds_x_left,2:3)];
Y_shift_rl = [coordinates(inds_y_righ,1) coordinates(inds_y_righ,2)-box_y coordinates(inds_y_righ,3)];
Y_shift_lr = [coordinates(inds_y_left,1) coordinates(inds_y_left,2)+box_y coordinates(inds_y_left,3)];
Z_shift_rl = [coordinates(inds_z_righ,1:2) coordinates(inds_z_righ,3)-box_z];
Z_shift_lr = [coordinates(inds_z_left,1:2) coordinates(inds_z_left,3)+box_z];

%%% corners
XY_shift_rtlb = [coordinates(inds_xy_rightop,1)-box_x coordinates(inds_xy_rightop,2)-box_y coordinates(inds_xy_rightop,3)];
XY_shift_rblt = [coordinates(inds_xy_righbot,1)-box_x coordinates(inds_xy_righbot,2)+box_y coordinates(inds_xy_righbot,3)];
XY_shift_ltrb = [coordinates(inds_xy_lefttop,1)+box_x coordinates(inds_xy_lefttop,2)-box_y coordinates(inds_xy_lefttop,3)];
XY_shift_lbrt = [coordinates(inds_xy_leftbot,1)+box_x coordinates(inds_xy_leftbot,2)+box_y coordinates(inds_xy_leftbot,3)];

XZ_shift_rtlb = [coordinates(inds_xz_rightop,1)-box_x coordinates(inds_xz_rightop,2) coordinates(inds_xz_rightop,3)-box_z];
XZ_shift_rblt = [coordinates(inds_xz_righbot,1)-box_x coordinates(inds_xz_righbot,2) coordinates(inds_xz_righbot,3)+box_z];
XZ_shift_ltrb = [coordinates(inds_xz_lefttop,1)+box_x coordinates(inds_xz_lefttop,2) coordinates(inds_xz_lefttop,3)-box_z];
XZ_shift_lbrt = [coordinates(inds_xz_leftbot,1)+box_x coordinates(inds_xz_leftbot,2) coordinates(inds_xz_leftbot,3)+box_z];

YZ_shift_rtlb = [coordinates(inds_yz_rightop,1) coordinates(inds_yz_rightop,2)-box_y coordinates(inds_yz_rightop,3)-box_z];
YZ_shift_rblt = [coordinates(inds_yz_righbot,1) coordinates(inds_yz_righbot,2)-box_y coordinates(inds_yz_righbot,3)+box_z];
YZ_shift_ltrb = [coordinates(inds_yz_lefttop,1) coordinates(inds_yz_lefttop,2)+box_y coordinates(inds_yz_lefttop,3)-box_z];
YZ_shift_lbrt = [coordinates(inds_yz_leftbot,1) coordinates(inds_yz_leftbot,2)+box_y coordinates(inds_yz_leftbot,3)+box_z];

%%%

%coord_ext = [coordinates; X_shift_rl; X_shift_lr;Y_shift_rl;Y_shift_lr;Z_shift_rl;Z_shift_lr]; % extended coordinate list
coord_ext = [coordinates; X_shift_rl; X_shift_lr;Y_shift_rl;Y_shift_lr;Z_shift_rl;Z_shift_lr;...
    XY_shift_rtlb; XY_shift_rblt; XY_shift_ltrb; XY_shift_lbrt; XZ_shift_rtlb; XZ_shift_rblt; XZ_shift_ltrb; ...
    XZ_shift_lbrt; YZ_shift_rtlb; YZ_shift_rblt; YZ_shift_ltrb; YZ_shift_lbrt]; % extended coordinate list

%% finding 4 closest neighbors and their distances with each atom (the closes is the atom itself with dis = 0) 
[neigh_list, neigh_dist] = knnsearch(coord_ext,coord_ext,'k',5);

neigh_list(neigh_dist>cutoff) = 0; % consideing neighbours closser than cutoff 
neigh_list = neigh_list(:,2:5); % first column (the closest atoms) are the atoms themselves

neigh_dir = neigh_dir_func(coord_ext,neigh_list); % finding the bond directions for each atom
neigh_num = sum(logical(neigh_list),2); % finding number of neighbours for each atom
end

function neigh_dir = neigh_dir_func(coord,neigh_list)

N = size(coord,1); % number of atoms
neigh_dir = zeros(3,4,N); % 3d tensor of vectors to 4 neighbours

for i = 1:N % for all atoms
    for j = 1:4 % for all 4 neighbours
        if neigh_list(i,j)~=0
            neigh_dir(:,j,i) = coord(neigh_list(i,j),:) - coord(i,:);
        end
    end
end

neigh_dir = neigh_dir./vecnorm(neigh_dir+eps,2,1);  % normalizing the norm of directions

end
