clear all;close all;clc;
restoredefaultpath;
addpath('/ImagePTE1/ajoshi/code_farm/svreg/src');
addpath('/ImagePTE1/ajoshi/code_farm/svreg/3rdParty/')
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

%%
r=readdfs('/ImagePTE1/ajoshi/code_farm/bfp/supp_data/bci32kright.dfs');
l=readdfs('/ImagePTE1/ajoshi/code_farm/bfp/supp_data/bci32kleft.dfs');

c=jet(10);
for j=1:10
    f=view_patch(r);
    view(90,0);
    %mysphere(cursor_info1.Position,5,'r',100)
    mysphere([109.2022   48.3260  156.7884],4,c(j,:),100);
    material dull;
    axis tight;
    saveas(f,sprintf('brain_%d.png',j));
    close all;    
end

%%
load('/ImagePTE1/ajoshi/code_farm/bfp/supp_data/USCLobes_grayordinate_labels.mat');
r=readdfs('/ImagePTE1/ajoshi/code_farm/bfp/supp_data/bci32kright_smooth.dfs');
c=0.5+0*r.vertices;
c(labels(32492+1:2*32492)==200,1:2)=1;
patch('vertices',r.vertices,'faces',r.faces,'facevertexcdata',c,'edgecolor','none','facecolor','flat');
view(90,0);axis off;axis tight;axis equal;
camlight;
material dull;

