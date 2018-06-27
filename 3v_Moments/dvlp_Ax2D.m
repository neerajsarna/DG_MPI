function Ax = dvlp_Ax2D(M)
all_idx = cell(M+1,1);

%% we collect all the even and odd basis indices
for i = 0 : M
    [all_idx{i+1},~,~] = IDX_Full(i);
end

total_var = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:end),'Un',0)));

ja = [];
va = [];

% the id of all the degrees apart from M-1
iaOld = 1:sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:end-1),'Un',0)));

% loop over all the lower order moments
for i = 0 : M-1
    id_temp = all_idx{i+1};
    id_temp(:,1) = id_temp(:,1) + 1;
    
    % iab is the id of the rows in all_index{i+2} which are the same as
    % id_temp
    
    [~,iab,~] = intersect(all_idx{i+2},id_temp,'rows');
    
    % the result if not necessarily in the correct order
    iab = sort(iab);
    
    % we need to add the shift to account for all the previous variables
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:i+1),'Un',0)));
    ja = [ja iab' + shift_id];
    va = [va sqrt(id_temp(:,1))'];
end

% account for the symmetricity of the matrix
ia = [iaOld, ja];
ja = [ja, iaOld];
va = [va, va];

Ax = sparse(ia,ja,va,total_var,total_var);
end
