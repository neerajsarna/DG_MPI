function B = dvlp_BInflow2D(M)

[id_Odd,id_Even] = get_id_Odd(M);

n_odd = sum(cell2mat(cellfun(@(a) length(a),id_Odd,'Un',0)));
n_even = sum(cell2mat(cellfun(@(a) length(a),id_Even,'Un',0)));
n_tot = n_odd + n_even;

all_idx = cell(M+1,1);

%% we collect all the even and odd basis indices
for i = 0 : M
    [all_idx{i+1},~,~] = IDX_Full(i);
end

% idx of the Even and the odd Functions
all_idx_Even = cellfun(@(a) a(rem(a(:,1),2)==0,:),all_idx,'Un',0);
all_idx_Odd = cellfun(@(a) a(rem(a(:,1),2)~=0,:),all_idx(1:end),'Un',0);

%% develop the matrices for 1D
% M is the maximum possible tensor degree
odd_1D = 1:2:M;
even_1D = 0:2:M;

Boe_1D = zeros(length(odd_1D),length(even_1D));

for i = 1:length(odd_1D)
    for j = 1:length(even_1D)
        Boe_1D(i,j) = 2 * HermiteHalfSpace(odd_1D(i),even_1D(j));
    end    
end

%% develop for 2D
ia = [];
ja = [];
va = [];

for i = 1 : length(all_idx_Odd)
    % sum of all the previous ones
    shift_ia = sum(cell2mat(cellfun(@(a) size(a,1),all_idx_Odd(1:i-1),'Un',0)));
    
   
    for j = 1 : length(all_idx_Even)
        
        shift_ja = sum(cell2mat(cellfun(@(a) size(a,1),all_idx(1:j-1),'Un',0)));
   
        for k = 1 : size(all_idx_Odd{i},1)
         y_Odd = all_idx_Odd{i}(k,2); 
         z_Odd = all_idx_Odd{i}(k,3);
         
            for l = 1 : size(all_idx_Even{j},1)
                % if the y and z components are the same, only then do the
                % computations
                if y_Odd == all_idx_Even{j}(l,2)
                   if z_Odd == all_idx_Even{j}(l,3)
                        ia = [ia shift_ia + k];
                        ja = [ja id_Even{j}(l)];
                        va = [va -Boe_1D(ceil(all_idx_Odd{i}(k,1)/2),ceil((all_idx_Even{j}(l,1) + 1)/2))];
                   end
                end
            end
        end
    end
end

B = sparse(ia,ja,va,n_odd,n_tot);

ia = flatten_cell(id_Odd);
B(:,ia) = eye(n_odd);
end