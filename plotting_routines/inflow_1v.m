clear all;
range = 0:1;

num_proc = 1;

M_values = [4,5,6,7,15];
theta_values = cell(length(M_values),1);
for i = 1 : length(M_values)
    [X,theta_values{i}] = get_moment(1,M_values(i),5);
end

plot(X,theta_values{1},X,theta_values{2},X,theta_values{3},...
     X,theta_values{4},X,theta_values{5});
legend('4','5','6','7','15','location','best');
xlim([0.75 1.0]);

function [X,quantity] = get_moment(num_proc,M,id)

X = [];
quantity = [];

for j = 0 : num_proc-1
filename = strcat('../2x1v_Inflow/M',num2str(M),...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

X = [X;field1(:,1)];

quantity = [quantity;field1(:,id)];

end

end



