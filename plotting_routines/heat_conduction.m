clear all;
range = 0:1;

num_proc = 8;

M = 8;
M2 = 5;
M3 = 4;
Nx = 320;

X = [];
Y = [];

XM3 = [];
YM3 = [];

theta = [];
rho = [];

thetaM2 = [];
rhoM2 = [];

thetaM3 = [];
rhoM3 = [];
    
for j = 0 : num_proc-1
filename = strcat('../2x3v_moments_HC/M',num2str(M),...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

X = [X;field1(:,1)];
Y = [Y;field1(:,2)];

IDtheta = [6,8,9];
theta = [theta;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3];

rho = [rho; field1(:,3)];

end

for j = 0 : num_proc-1
filename = strcat('../2x3v_moments_HC/M',num2str(M2),...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

IDtheta = [6,8,9];
thetaM2 = [thetaM2;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3];

rhoM2 = [rhoM2; field1(:,3)];
end

for j = 0 : num_proc-1
filename = strcat('../2x3v_moments_HC/M',num2str(M3),...
                  '/result',num2str(j),'_Kn_0p1.txt');    
             
field1 = dlmread(filename,'\t');

XM3 = [XM3;field1(:,1)];
YM3 = [YM3;field1(:,2)];

IDtheta = [6,8,9];
thetaM3 = [thetaM3;sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        + field1(:,IDtheta(3)))/3];

rhoM3 = [rhoM3; field1(:,3)];
end

foldername= '../../StaRMAP_ver1p9/SBP_Traditional/generic_1D/result_HC2D_1x3v_Kn0p1_theta1/hc_tend_1_points_300_neqn_';
filename = strcat(foldername,num2str(45),'.txt');    
field = dlmread(filename,'\t');

X2 = field(1,:);
IDtheta = [4,5,6];                    
thetaM45 = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
                        +field(IDtheta(3),:))/3;
rhoM45 = field(2,:);
     
filename = strcat(foldername,num2str(M),'.txt');    
field = dlmread(filename,'\t');
theta_sbpM = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
                        +field(IDtheta(3),:))/3;
rho_sbpM = field(2,:);                  

filename = strcat(foldername,num2str(M3),'.txt');    
field = dlmread(filename,'\t');
theta_sbpM3 = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
                        +field(IDtheta(3),:))/3;
                    
rho_sbpM3 = field(2,:);                  

[x,theta] = compute_theta(X,Y,theta,Nx);
[~,thetaM2] = compute_theta(X,Y,thetaM2,Nx);
[xM3,thetaM3] = compute_theta(XM3,YM3,thetaM3,500);

[~,rho] = compute_theta(X,Y,rho,Nx);
[~,rhoM2] = compute_theta(X,Y,rhoM2,Nx);
[~,rhoM3] = compute_theta(XM3,YM3,rhoM3,500);

%plot(X,theta,'-',X2,theta2,'-',X2,theta3,'-');
% plot(x,theta,'-',...
%      X2,theta_sbpM,'-',...
%      xM3,thetaM3,'-',...
%      X2,theta_sbpM3,'-',...
%      X2,thetaM45,'-');
 
plot(x,rho,'-',...
     X2,rho_sbpM,'-',...
     xM3,rhoM3,'-',...
     X2,rho_sbpM3,'-',...
     X2,rhoM45,'-');
 
legend('M=8','M=8sbp','M=4',...
       'M=4sbp','M45','location','best');
%xlim([0.5 1]);

% develop x y and theta from the 2D data
function [x,theta] = compute_theta(X,Y,theta,Nx)
[~,id] = sort(Y);
X = X(id);
x = X(1:Nx);
theta = theta(id);
theta = theta(1:Nx);
end



