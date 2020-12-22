function [S, E] = temp_matching_func(neigh_dirs,neigh_list,neigh_num,atoms_num,score_min)
% calculating the similariy score between structure of water molecules in the box and the
% staggered and eclipced templates

Pp = perms([1 2 3]);
I = eye(3);
P = zeros(6,9);
for i = 1:6
    It = I(:,Pp(i,:));
    P(i,:) = It(:)';
end

Ts = [-1 0.5 0.5; 0.5 -1 0.5; 0.5 0.5 -1]; % staggered template
Te = [0.5 -0.5 -0.5; -0.5 0.5 -0.5; -0.5 -0.5 0.5]; % eclipsed template

E = zeros(atoms_num,4);
S = zeros(atoms_num,4);

for i = 1:atoms_num
    for j = 1:4
        jth_neigh_i =neigh_list(i,j);
        if jth_neigh_i>i
            
            [Tx,discon_flag] = UtV_func(neigh_dirs(:,:,i),neigh_dirs(:,:,jth_neigh_i));
            nU = neigh_num(i)-1;
            nV = neigh_num(jth_neigh_i)-1;
            
            if  discon_flag==0
                [s, e] = score_func(Tx,Ts,Te,nU,nV,P);
            else
                neigh_list(i,j)=0;
                neigh_num(i) = neigh_num(i) -1;
                s = 0;
                e = 0;
            end
            
            S(i,j) = s;
            E(i,j) = e;
            if jth_neigh_i<=atoms_num
                S(jth_neigh_i,neigh_list(jth_neigh_i,:)==i) = s;
                E(jth_neigh_i,neigh_list(jth_neigh_i,:)==i) = e;
            end
            
        end
    end
end

E(E<S)=0;
S(S<E)=0;
E(E<score_min)=0;
S(S<score_min)=0;

end

function [Tx,disconnection_flag] = UtV_func(U,V)
% calculating U'*V
Tx = U'*V;
[ii, jj] = find(abs(Tx+1)<1e-7);
if isempty(ii) % if there is no central O-O bond
    disconnection_flag = 1; % there is no tetrahedral structure (disconnection in O-O)
else
    Tx(ii(1),:)=[];
    Tx(:,jj(1))=[];
    disconnection_flag = 0;
end
end


function [s, e] = score_func(Tx,Ts,Te,nU,nV,P)
DS = zeros(3,3);
DE = DS;
lamb=0.15;
if nV>0 && nU>0
    DS(1:nV,1:3) = dist(Tx(1:nU,1:nV)',Ts(1:nU,:)).^2;
    DE(1:nV,1:3) = dist(Tx(1:nU,1:nV)',Te(1:nU,:)).^2;
    
    ds = min(P*DS(:))/(lamb*nV*nU);
    de = min(P*DE(:))/(lamb*nV*nU);
    
    s = exp(-ds);
    e = exp(-de);
else
    s=0;
    e=0;
end
end

