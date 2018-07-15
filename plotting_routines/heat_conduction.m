clear all;
range = 0:1;

num_proc = 8;

M = 3;
M2 = 4;
Nx = 300;

X = [];
Y = [];

theta = [];
rho = [];

thetaM2 = [];
rhoM2 = [];


for j = 0 : 0
filename = strcat('../2x3v_moments_HC_Adp/M_',num2str(M),...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

X = [X;field1(:,1)];
Y = [Y;field1(:,2)];

IDtheta = [6,8,9];
theta = [theta;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3];

rho = [rho; field1(:,3)];

end

for j = 0 : 0
filename = strcat('../2x3v_moments_HC_Adp/M_',num2str(M2),...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

IDtheta = [6,8,9];
thetaM2 = [thetaM2;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3];

rhoM2 = [rhoM2; field1(:,3)];
end

%% sbp results
foldername= '../../StaRMAP_ver1p9/SBP_Traditional/generic_1D/result_HC2D_1x3v_Kn0p1_theta1/hc_tend_1_points_300_neqn_';     
IDtheta = [4,5,6];
filename = strcat(foldername,num2str(M),'.txt');    
field = dlmread(filename,'\t');
theta_sbpM = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
                        +field(IDtheta(3),:))/3;
rho_sbpM = field(2,:);                  

filename = strcat(foldername,num2str(M2),'.txt');    
field = dlmread(filename,'\t');
theta_sbpM2 = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
                        +field(IDtheta(3),:))/3;
                    
rho_sbpM2 = field(2,:);                  


plot(field(1,:),theta_sbpM,...
     field(1,:),theta_sbpM2,'-');

legend('M3sbp','M4sbp');
xlim([0.85 1]);
 

% develop x y and theta from the 2D data
function [x,theta] = compute_theta(X,Y,theta,Nx)
[~,id] = sort(Y);
X = X(id);
x = X(1:Nx);
theta = theta(id);
theta = theta(1:Nx);
end



