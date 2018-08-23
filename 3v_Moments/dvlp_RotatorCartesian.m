function Rotator = dvlp_RotatorCartesian(M,to_check)

num_var = 0;

% by idx_nt we denote the multi-index for the coordinate system with normal
% along n and tangent along t. 

% x positive and y positive
% bc_id = 1 and bc_id 3
all_idx_x_y = cell(M+1,1);
all_idx_y_x = cell(M+1,1);

%% we collect all the even and odd basis indices
for i = 0 : M
    [all_idx_x_y{i+1},~,~] = IDX_Full(i);
    
    % swap x and y, z remain the same
    all_idx_y_x{i+1}(:,1) = all_idx_x_y{i+1}(:,2);
    all_idx_y_x{i+1}(:,2) = all_idx_x_y{i+1}(:,1);
    all_idx_y_x{i+1}(:,3) = all_idx_x_y{i+1}(:,3);

    
    num_var = num_var + size(IDX_Full(i),1);
end

%% allocate memory for rotator
% four different boundaries
Rotator = cell(4,1);

for i = 1 : 4
    Rotator{i} = zeros(num_var,num_var);
end

%% bc_id = 1, x = 1
Rotator{1} = eye(num_var);

%% bc_id = 3, x = 0 
x_indices = cell2mat(cellfun(@(a) a(:,1),all_idx_x_y(:),'Un',0));
y_indices = cell2mat(cellfun(@(a) a(:,2),all_idx_x_y(:),'Un',0));

for i = 1 : size(x_indices)
    Rotator{3}(i,i) = (-1)^(x_indices(i)) * (-1)^(y_indices(i));
end

%% bc_id = 2 and bc_id = 4
xyz_indices = cell2mat(cellfun(@(a) a(:,[1 2 3]),all_idx_x_y(:),'Un',0)); 
yxz_indices = cell2mat(cellfun(@(a) a(:,[1 2 3]),all_idx_y_x(:),'Un',0)); 

% where does yx indices lie in the xy indices 
for i = 1:length(yxz_indices)
    Lia = ismember(xyz_indices,yxz_indices(i,:),'rows');
    % negative in the x direction
    Rotator{2}(i,Lia) = (-1)^(xyz_indices(Lia,1));
    % negative in the y direction
    Rotator{4}(i,Lia) = (-1)^(xyz_indices(Lia,2));
end

% choose to check the rotator
if to_check
    check_rotator(Rotator,xyz_indices,yxz_indices);
end
end

function [] = check_rotator(Rotator,xyz,yxz)

norm1 = norm(Rotator{1}*xyz-xyz);
norm2 = norm(abs(Rotator{2}*xyz)-yxz);
norm3 = norm(abs(Rotator{3}*xyz)-xyz);
norm4 = norm(abs(Rotator{4}*xyz)-yxz);


if norm1 ~= 0 || norm2 ~= 0 || norm3 ~= 0 || norm4 ~= 0
    error('incorrect rotator');
end

end