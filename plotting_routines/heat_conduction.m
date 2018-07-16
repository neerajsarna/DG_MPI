clear all;
range = 0:1;

num_proc = 8;

%% moment solution

M = 3;
M2 = 5;
Nx = 300;

X = [];
Y = [];

theta = [];
rho = [];
sigma11 = [];
qx = [];

thetaM2 = [];
rhoM2 = [];

thetaAdp = [];
rhoAdp = [];


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

for j = 0 : 0
filename = strcat('../2x3v_moments_HC_Adp/M_3_7',...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

IDtheta = [6,8,9];
thetaAdp = [thetaAdp;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3];

rhoAdp = [rhoAdp; field1(:,3)];
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


plot(X,(theta-thetaM2));

disp('error M3');
disp(sqrt(trapz(X,abs(theta-thetaM2).^2)));

disp('error M37');
disp(sqrt(trapz(X,abs(thetaAdp-thetaM2).^2)));

legend('error M3','error M37');
xlim([0 1]);
 

% develop x y and theta from the 2D data
function [x,theta] = compute_theta(X,Y,theta,Nx)
[~,id] = sort(Y);
X = X(id);
x = X(1:Nx);
theta = theta(id);
theta = theta(1:Nx);
end



