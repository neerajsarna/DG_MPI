function [] = develop_System(M)

disp('value of M: ');
disp(M);

Ax = dvlp_Ax1D(M);
P = -eye(M);
P(1,1) = 0;
P(2,2) = 0;
P(3,3) = 0;

Binflow = dvlp_BInflow1D(M);
Bwall = dvlp_BWall1D(M);

Binflow =  stabilize_boundary(Ax,Binflow);
Bwall = stabilize_boundary(Ax,Bwall);

penaltyInflow = dvlp_penalty(Ax,Binflow);
penaltyWall = dvlp_penalty(Ax,Bwall);

rotator = cell(4,1);

for i = 1 : 4
    rotator{i} = eye(M);
end

for i = 1 : M
    if mod(i,2) == 0
        rotator{3}(i,i) = -1;
    end
end


%% write Ax
data = get_sparse_data(Ax);
filename = strcat('Ax/Ax',num2str(M),'.txt');
mkdir Ax;
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write P
data = get_sparse_data(P);
mkdir P;
filename = strcat('P/P',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write Binflow
data = get_sparse_data(Binflow);
mkdir Binflow;
filename = strcat('Binflow/Binflow',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write Bwall
data = get_sparse_data(Bwall);
mkdir Bwall;
filename = strcat('Bwall/Bwall',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write penaltyInflow
data = get_sparse_data(penaltyInflow);
filename = strcat('Binflow/penalty_inflow',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write penaltyWall
data = get_sparse_data(penaltyWall);
filename = strcat('Bwall/penalty_wall',num2str(M),'.txt');
dlmwrite(filename,size(data,1),'precision',16);
dlmwrite(filename,data,'delimiter',' ','-append','precision',16);

%% write the rotator
mkdir Rotator;
for i = 1 : 4
    data = get_sparse_data(rotator{i});
    filename = strcat('Rotator/rotator',num2str(M),'_',num2str(i),'.txt');
    dlmwrite(filename,size(data,1),'precision',16);
    dlmwrite(filename,data,'delimiter',' ','-append','precision',16);
end

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


function B_ID1 = dvlp_B_ID1(varargin)

B_ID2 = varargin{1};

if nargin == 1
    % number of columns are the number of variables
    id_odd = 2:2:size(B_ID2,2);
end

if nargin == 2
    [id_odd,~] = get_id_Odd(varargin{2});
    id_odd = flatten_cell(id_odd);
end

    B_ID1 = B_ID2;
    
    % flip the signs of the odd variables
    B_ID1(:,id_odd) = -B_ID1(:,id_odd);
    
end


