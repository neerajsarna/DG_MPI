function P = dvlp_Prod2D(M)

if M < 3 
    error('M should be greater than 3');
end

n_eqn = 0;

for i = 0 : M
    [idx,~,~] = IDX_Full(i);
    n_eqn = n_eqn + size(idx,1);
end

% 3 for conservation laws and 6 for cross coupling
nnz = n_eqn - 3 + 6;

ia = zeros(1,nnz);
ja = zeros(1,nnz);
va = zeros(1,nnz);

ia(1:10) = [4,4,4,5,6,6,6,7,7,7];
ja(1:10) = [4,6,7,5,4,6,7,4,6,7];
va(1:10) = -[2/3,-1/3,-1/3, ...
            1, ...
            -1/3,2/3,-1/3, ...
            -1/3,-1/3,2/3];
        
% add the diagonal elements, minus sign to account for the relaxation
ia(11:nnz) =8:n_eqn;
ja(11:nnz) =8:n_eqn;
va(11:nnz) = -ones(1,length(8:n_eqn));

P = sparse(ia,ja,va,n_eqn,n_eqn);
end

