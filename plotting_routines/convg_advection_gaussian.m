clear all;

filename = '../2D_advection_gaussian/error.txt';
error = dlmread(filename,'\t');

order = -1;
[reference_order] = exact_order(sqrt(error(end,1)),error(end,2), ...
                                sqrt(error(1,1)),order);
                            
loglog(sqrt(error(:,1)),error(:,2),'-*', ...
          reference_order(1:2),reference_order(3:4),'-o');
legend('numerical','expected');
grid on;
xlim([min(sqrt(error(:,1))) max(sqrt(error(:,1)))]);