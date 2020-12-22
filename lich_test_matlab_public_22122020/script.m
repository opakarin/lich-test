% Author: Golnaz Roudsari
% Email: golnaz.roudsari@helsinki.fi

%%
 clear
% clc
%%
% filename = 'plannar_full_ice';
%   filename = 'npt';
%filename = 'cubic';

    filetype = '.xyz';

tic
[watoms_num, watoms_coord, satoms_num, satoms_coord, box_x, box_y, box_z] = read_file_func([filename filetype]);
 
%%
cutoff = 3.5; % a distance cutoff of the first nearest neighbour shell
[neigh_list, neigh_dir, neigh_num, inds_ext] = neigh_info_func(watoms_coord, box_x, box_y, box_z, cutoff);

%%
min_score = 0.5; % the least cosiderable similarity score between the structure of atoms in the box and a template 
[Stg_score, Ecl_score] = temp_matching_func(neigh_dir,neigh_list,neigh_num,watoms_num,min_score);

%%
Labels = ice_labelling_func(logical(Stg_score), logical(Ecl_score), neigh_list, neigh_num, inds_ext, watoms_num);
toc
%%
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
        fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','I',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
    else
      fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','L',watoms_coord(i,1),watoms_coord(i,2),watoms_coord(i,3));
    end
end
for i=1:satoms_num
       fprintf(fileID,'\n%s %12.8f %12.8f %12.8f ','X',satoms_coord(i,1),satoms_coord(i,2),satoms_coord(i,3));
end
fclose(fileID);
%%
disp( [ 'from ' num2str(watoms_num) ' oxygens:'  num2str(sum(Labels==1)) ' C; ' num2str(sum(Labels==2)) ' H; ' ...
    num2str(sum(Labels==3)) ' MI; ' num2str(sum(Labels==4)) ' CI; ' num2str(sum(Labels==5))...
    ' HI; ' num2str(sum(Labels==6)) ' I; '  num2str(sum(Labels==0)) ' L; '])

end
