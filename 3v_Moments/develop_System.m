function [] = develop_System(M)

disp('value of M: ');
disp(M);

Ax = dvlp_Ax2D(M);
P = dvlp_Prod2D(M);
Binflow = dvlp_BInflow2D(M);
Bwall = dvlp_BWall2D(M);

Binflow =  stabilize_boundary(Ax,Binflow,M);
Bwall =  stabilize_boundary(Ax,Bwall,M);

rotator = dvlp_RotatorCartesian(M,false);

penaltyInflow = dvlp_penalty_char(Ax,Binflow);
penaltyWall = dvlp_penalty_char(Ax,Bwall);


%% write Ax
data = get_sparse_data(Ax);
filename = strcat('Ax/Ax',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write P
data = get_sparse_data(P);
filename = strcat('P/P',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write Binflow
data = get_sparse_data(Binflow);
filename = strcat('Binflow/Binflow',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write penaltyInflow
data = get_sparse_data(penaltyInflow);
filename = strcat('Binflow/penalty_inflow',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write Bwall
data = get_sparse_data(Bwall);
filename = strcat('Bwall/Bwall',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write penaltyWall
data = get_sparse_data(penaltyWall);
filename = strcat('Bwall/penalty_wall',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write the rotator
for i = 1 : 4
data = get_sparse_data(rotator{i});
filename = strcat('Rotator/rotator',num2str(M),'_',num2str(i),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);
end
end

function n_eqn = compute_size(M)
n_eqn = (3 * M^2 -3*M+8)/2;
end
function [data] = get_sparse_data(mat)
[ii,jj,va] = find(mat);

% change to c format
ii = ii - 1;
jj = jj - 1;

data = zeros(length(ii),3);
data(:,1) = ii;
data(:,2) = jj;
data(:,3) = va;
end

function [penalty] = dvlp_penalty_char(An,B)
modAn = compute_Amod(An);
Xminus = compute_Xminus(An);

if norm(full(B)) == 0
    penalty = Xminus * 0;
else
    penalty = 0.5 * (An-modAn) * Xminus/(B*Xminus);
end
end

function modA = compute_Amod(A)

% eig does not support sparse matrices
[V,D] = eig(full(A));
modA = V * abs(D)/V;
end

% compute the eigenvector corresponding to the negative eigenvalues
function Xminus = compute_Xminus(A)
[V,D] = eig(full(A));
D= D(sub2ind(size(D),1:size(D,1),1:size(D,2)));
loc_neg = D<-(1e-10);
Xminus = V(:,loc_neg);
end


