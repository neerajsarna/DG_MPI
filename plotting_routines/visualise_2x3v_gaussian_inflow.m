clear all;

filename = '../2x3v_moments_gaussian_inflow/result0.txt';
field = dlmread(filename,'\t');

X = field(:,1);
Y = field(:,2);
Z = field(:,3);

figure(3);
scatter3(X,Y,Z);