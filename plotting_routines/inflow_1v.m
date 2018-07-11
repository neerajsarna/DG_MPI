clear all;
range = 0:1;

num_proc = 1;
M_values = [4:1:15];

theta_values = cell(length(M_values),1);
for i = 1 : length(M_values)
    if M_values(i) >= 11 && M_values(i) <= 14
        num_proc = 2;
    else 
        num_proc = 1;
    end
    
        [X,theta_values{i}] = get_moment(num_proc,M_values(i),5);
end

error = [];
% for i = 1 : (length(M_values))
%     plot(X,theta_values{i},'LineWidth',2);
%     xlim([0.85 1.0]);
%     hold on;
% end

for i = 1 : (length(M_values)-1)
    error = [error, compute_error(X,theta_values{i},theta_values{end})];
end

loglog(M_values(1:end-1),error);
xlim([4 15]);


% [~, hobj, ~, ~] = legend('5','7','9','15','location','best');
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',2)

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

function error = compute_error(X,data1,data2)
error = (data1-data2).*(data1-data2);
error = sqrt(trapz(X,error));
end



