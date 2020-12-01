clear all;close all;clc;
restoredefaultpath;
addpath('/ImagePTE1/ajoshi/code_farm/svreg/src');

figure;
mysphere([0,0,0],100,[.5,.5,.5],100);
hold on;axis equal;material dull;camlight;


cml=jet(125);


r =100;
phi=.3*rand(125,1);
theta=2*pi*rand(125,1);
%[phi,theta]=meshgrid(phirange,thetarange);
x = r .* sin(phi(:)) .* cos(theta(:));
y = r .* cos(phi(:));
z = r .* sin(phi(:)) .* sin(theta(:));

for j=1:length(x)
    mysphere([x(j,:),y(j,:),z(j,:)],3,cml(j,:),10);
end

view(180,-60);


r=120*(1+phi).^3;
x = r .* sin(phi(:)) .* cos(theta(:));
y = (r) .* cos(phi(:));
z = r .* sin(phi(:)) .* sin(theta(:));

for j=1:length(x)
    mysphere([x(j,:),y(j,:),z(j,:)],3,cml(j,:),10);
end
camlight;




