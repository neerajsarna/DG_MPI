function B = dvlp_BWall1D(M)

BInflow = dvlp_BInflow1D(M);

eqW = BInflow(1,:)/BInflow(1,1);

B = BInflow - BInflow(:,1)*eqW;

% first equation
B(:,2) = 0;
B(1,2) = 1;
end