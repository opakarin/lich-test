function [watoms_num, watoms_coord, satoms_num, satoms_coord, box_x,box_y,box_z] = read_file_func(filename)
%%%%% [num_waters,coordinate_waters,num_surface,coordinate_surface,box_size1,box_size2,box_size3]...
%%%%% ...= read_file_func(filename)

fileID = fopen(filename);
C = textscan(fgetl(fileID),'%d');
num_atoms = C{1};

box_size = textscan(fgetl(fileID),'%s %f %f %f');
box_x = box_size{2};
box_y = box_size{3};
box_z = box_size{4};


atoms_xyz = textscan(fileID, '%s %f %f %f %d %d %f %f %f');
Type = atoms_xyz{1};
X = atoms_xyz{2};
Y = atoms_xyz{3};
Z = atoms_xyz{4};

winds = strcmp(Type,'O');
sinds = strcmp(Type,'S');
watoms_coord = [X(winds) Y(winds) Z(winds)];
satoms_coord = [X(sinds) Y(sinds) Z(sinds)];
watoms_num = size(watoms_coord,1);
satoms_num = size(satoms_coord,1);

fclose(fileID);
end


