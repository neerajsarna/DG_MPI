function f = dvlp_Ax1D(n_eqn)
var = 0:1:n_eqn-1;

id_odd = 2:2:n_eqn;
id_even = 1:2:n_eqn;
n_odd = length(id_odd);
n_even = length(id_even);

Aoe = zeros(n_odd,n_even);
Aoe(1:n_odd+1:end) = sqrt(var(id_odd));
Aoe((n_odd+1):n_odd+1:end) = sqrt(var(id_odd(1:(end-mod(n_eqn+1,2))))+1);

A = zeros(n_eqn,n_eqn);

A(id_odd,id_even) = Aoe;
A(id_even,id_odd) = Aoe';

f = A;
end