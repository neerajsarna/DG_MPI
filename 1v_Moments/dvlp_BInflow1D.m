function f = dvlp_BInflow1D(n_eqn)

var = 0:1:n_eqn-1;

id_odd = 2:2:n_eqn;
id_even = 1:2:n_eqn;

n_odd = length(id_odd);
n_even = length(id_even);

Boe = zeros(n_odd,n_even);

for i = 1 : n_odd
    for j = 1 : n_even
        Boe(i,j) = 2 * HermiteHalfSpace(var(id_odd(i)),var(id_even(j)));
    end
end

B = zeros(length(id_odd),length(id_even));

B(:,id_odd) = eye(length(id_odd));
B(:,id_even) = -Boe;

f = B;
end