clear all;

filename = '../1D_advection_sin/result0.txt';
field = dlmread(filename,'\t');

X = field(:,1);
Y = field(:,2);
Z = field(:,3);

t = 1.0;
%Zexact = exp(-(X-t-0.5).^2*100);
Zexact = cos(pi * (X-t));

figure(1);
scatter3(X,Y,abs(Z-Zexact));

figure(3);
scatter3(X,Y,Z);

filename = '../1D_advection_sin/error_sin.txt';
error = dlmread(filename,'\t');

% we do not take the square root because we only change the cells in the 
% x-direction
order = -1;
[reference_order] = exact_order(error(end-1,1),error(end-1,2), ...
                                error(1,1),order);
          
figure(2);
loglog(error(:,1),error(:,2),'-*', ...
          reference_order(1:2),reference_order(3:4),'-o');
legend('numerical','expected');
grid on;
xlim([min(error(:,1)) max(error(:,1))]);