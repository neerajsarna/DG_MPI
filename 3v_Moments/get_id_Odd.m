% return the global ids of odd variables
function [id_Odd,id_Even] = get_id_Odd(M)
all_idx = cell(M+1,1);

%% we collect all the even and odd basis indices
% the triplet for every different moment order
for i = 0 : M
    [all_idx{i+1},~,~] = IDX_Full(i);
end

%% ids of the odd variables
id_Odd = cell(M+1,1);
id_Even = cell(M+1,1);

for i = 0: M
    
    shift_id = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:i),'Un',0)));
    
    id_temp = all_idx{i+1}(rem(all_idx{i+1}(:,1),2) ~= 0,:);
    [~,ia,~] = intersect(all_idx{i+1},id_temp,'rows');
    
    ia = sort(ia);
    % we directly store the global index along with the shift. 
    id_Odd{i+1} = ia' + shift_id;
    
    id_temp = all_idx{i+1}(rem(all_idx{i+1}(:,1),2) == 0,:);
    [~,ia,~] = intersect(all_idx{i+1},id_temp,'rows');
    
    ia = sort(ia);
    % we store the global index along with the shift
    id_Even{i+1} = ia' + shift_id;
end
end