function Labels = ice_labelling_func(Stg_score, Ecl_score, neigh_list, neigh_num, inds_ext, watoms_num)

Labels = zeros(watoms_num,1);
Labels(sum(Stg_score,2)==4 & sum(Ecl_score,2)==0,1)=1; % C
Labels(sum(Stg_score,2)==3 & sum(Ecl_score,2)==1 & sum(Ecl_score+Stg_score,2)==4,1)=2; % H

% list = zeros(watoms_num,1);
% list(sum(Stg_score,2)+sum(Ecl_score,2)>0 & Labels==0,1)=1; % I
% 
% for i = find(list)'
%    i_neighs = neigh_list(i,1:neigh_num(i));
%    i_neighs(i_neighs>watoms_num) = inds_ext(i_neighs(i_neighs>watoms_num)-watoms_num);
%    
%    if sum(Stg_score(i,:).*Ecl_score(i,:))>0
%         if sum(Labels(i_neighs,1)==1)==3
%            Labels(i,1) = 1; % C
%         elseif sum(Labels(i_neighs,1)==2)==3
%            Labels(i,1) = 2; % H 
%         end
%    end              
% end

list = zeros(watoms_num,1);
list(sum(Stg_score+Ecl_score,2)>0 & Labels==0,1)=1; % I

for i = find(list)'
   i_neighs = neigh_list(i,1:neigh_num(i));
   i_neighs(i_neighs>watoms_num) = inds_ext(i_neighs(i_neighs>watoms_num)-watoms_num);
   
   if any(Labels(i_neighs,1)==1) && any(Labels(i_neighs,1)==2)
       Labels(i,1) = 3; %Mixed
   elseif any(Labels(i_neighs,1)==1)
       Labels(i,1) = 4; %IC
   elseif any(Labels(i_neighs,1)==2)
       Labels(i,1) = 5; %IH       
   end   
end

Labels(sum(Ecl_score,2)>=3,1)=6; % CH

list = zeros(watoms_num,1);
list(Labels==0,1)=1;
for i = find(list)'
   i_neighs = neigh_list(i,1:neigh_num(i));
   i_neighs(i_neighs>watoms_num) = inds_ext(i_neighs(i_neighs>watoms_num)-watoms_num);
   
   if any( (Labels(i_neighs,1)==3 | Labels(i_neighs,1)==4 | Labels(i_neighs,1)==5 | Labels(i_neighs,1)==6)' .* ( Ecl_score(i,1:neigh_num(i))+Stg_score(i,1:neigh_num(i)) ) )
       Labels(i,1) = 7; %I 
   end          
end

