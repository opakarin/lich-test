function Res = Main_Ice_Recog_func(filename,filetype,varargin)
%% This algorithm identifies different types of ice crystals in an MD simulation box
%
% Inputs:
%   1- filename (without extension)
%   2- filetype (extension, e.g., .xyz)
%   3- opts: optional parameters
%                - opts.cutoff: a distance cutoff of the first
%                  nearest neighbour shell (default: 3.5)
%                - opts.min_score: the least considerable similarity
%                  score between the structure of atoms in the box and
%                  a template (default: 0.6)
%
%
%    The first line in the input file should be the number of atoms. 
%    The second line should be one word (e.g., the name of the
%    file) followed by the sizes of the simiulation box in x, y and z 
%    directions, respectively, all in Angstroms. 
%    An example of a input file is given as follows.
%
%                     12000
%                     box_name 12.05 13.1 14.00    
%                     <element> <x> <y> <z>
%                     <element> <x> <y> <z>
%                     ...
%
%     In the input file, the capital letter 'O' is used for oxygen atoms,
%     and the rest of atoms are represented by the capital letter 'S'.
%
%
% Outputs:
%   1- The labeled file with the name filename_labeled.filetype
%   2- Res.runtime: the execusion time (only ice recognition algorithm)
%   3- Res.Ice_Types_nums: number of oxygen atoms with each ice structure
%           Ice structures:
%                           - C     (cubic)
%                           - H     (hexagonal)
%                           - CH    (clathrate-hydrates)
%                           - CI    (cubic-interfacial)
%                           - HI    (hexagonal-interfacial)
%                           - MI    (mixed-interfacial)
%                           - I     (interfacial)
%                           - L     (liquid)
%
%
% Author: Golnaz Roudsari
% Email: golnaz.roudsari@helsinki.fi

%% Parameters

if nargin<3
    opts = [];
end
if isfield(opts,'cutoff')
    cutoff = opts.cutoff;
else
    cutoff = 3.5;
end

if isfield(opts,'min_score')
    min_score = opts.min_score;
else
    min_score = 0.5;
end

addpath('./codes')

%% Main steps
disp('Ice structures identification ... ')
t = tic;
[watoms_num, watoms_coord, satoms_num, satoms_coord, box_x, box_y, box_z] = read_file_func([filename filetype]);
[neigh_list, neigh_dir, neigh_num, inds_ext] = neigh_info_func(watoms_coord, box_x, box_y, box_z, cutoff);
[Stg_score, Ecl_score] = temp_matching_func(neigh_dir,neigh_list,neigh_num,watoms_num,min_score);
Labels = ice_labelling_func(logical(Stg_score), logical(Ecl_score), neigh_list, neigh_num, inds_ext, watoms_num);
Res.runtime = toc(t);
disp('Identifican completed.')
%% Printing labeled file
    disp('Printing the labeled file ... ')
    fileID = fopen([filename '_labled' filetype],'w');
    fprintf(fileID,'%d',watoms_num+satoms_num);
    fprintf(fileID,'\n t=0');
    for i=1:watoms_num
        if Labels(i)==1
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','C',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==2
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','H',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==3
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','IM',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==4
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','IC',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==5
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','IH',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==6
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','CH',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==7
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','I',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        elseif Labels(i)==8
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','ICH',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        else
            fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','L',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
        end
    end
    for i=1:satoms_num
        fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','X',satoms_coord(i,1),satoms_coord(i,2),satoms_coord(i,3));
    end
    fclose(fileID);

%%
struct_type = {'Total'; 'C'; 'H'; 'CH'; 'CI'; 'HI'; 'CHI'; 'MI'; 'I'; 'L'};
struct_number = [watoms_num; sum(Labels==1); sum(Labels==2); sum(Labels==6); sum(Labels==4); sum(Labels==5); sum(Labels==8); sum(Labels==3); sum(Labels==7); sum(Labels==0)];
Res.Ice_Types_nums = table(struct_type, struct_number);

