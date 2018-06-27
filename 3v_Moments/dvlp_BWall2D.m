function B = dvlp_BWall2D(M)

BInflow = dvlp_BInflow2D(M);

eqW = BInflow(1,:)/BInflow(1,1);

B = BInflow - BInflow(:,1)*eqW;

% first equation
B(:,2) = 0;
B(1,2) = 1;
end