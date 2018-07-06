clear all;
range = 0:100:4700;

for i = 1 : length(range)
filename = strcat('../2x3v_moments_HC/result',num2str(range(i)),'.txt');    
field1 = dlmread(filename,'\t');


X = field1(:,1);
Y = field1(:,2);

IDtheta = [6,8,9];
theta = sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
                        +field1(:,IDtheta(3)))/3;

filename = '../../StaRMAP_ver1p9/SBP_Traditional/generic_1D/result_HC2D_1x3v_Kn0p1_theta1/hc_tend_1_points_300_neqn_3.txt';    
field = dlmread(filename,'\t');

X2 = field(1,:);
IDtheta = [4,5,6];                    
theta2 = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
                        +field(IDtheta(3),:))/3;
                    
figure(1);
plot(X,theta,'-o',X2,theta2,'-x');
drawnow;
end

% filename = '../2x3v_moments_HC/result0.txt';    
% field1 = dlmread(filename,'\t');
% 
% X = field1(:,1);
% Y = field1(:,2);
% 
% IDtheta = [6,8,9];
% theta = sqrt(2) * (field1(:,IDtheta(1))+field1(:,IDtheta(2))...
%                         +field1(:,IDtheta(3)))/3;
% 
% filename = '../../StaRMAP_ver1p9/SBP_Traditional/generic_1D/result_HC2D_1x3v_Kn0p1_theta1/hc_tend_1_points_300_neqn_3.txt';    
% field = dlmread(filename,'\t');
% 
% X2 = field(1,:);
% IDtheta = [4,5,6];                    
% theta2 = sqrt(2) * (field(IDtheta(1),:)+field(IDtheta(2),:)...
%                         +field(IDtheta(3),:))/3;
%                     
% plot(X,theta,'-o',X2,theta2,'-x');