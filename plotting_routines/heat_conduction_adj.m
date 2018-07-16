clear all;

%% adjoint moment solution

X = [];
XM = [];

theta = [];
vx = [];


for j = 0 : 0
filename = strcat('../2x3v_moments_HC_Adp/M_4Adj',...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

X = [X;field1(:,1)];

IDtheta = [6,8,9];
theta = [theta;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3];
end

for j = 0 : 0
filename = strcat('../2x3v_moments_HC_Adp/M_4',...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

XM = [XM;field1(:,1)];
vx = [vx;field1(:,4)];
end

filename = strcat('../2x3v_moments_HC_Adp/error_M3.txt');                 
error = dlmread(filename,'\t');

plot(X,abs(error(:,3)));
    
% develop x y and theta from the 2D data
function [x,theta] = compute_theta(X,Y,theta,Nx)
[~,id] = sort(Y);
X = X(id);
x = X(1:Nx);
theta = theta(id);
theta = theta(1:Nx);
end



