%Conformal Array Code (1D)
%Assumes the curvature is a segment of a circle centered at (0,0)
%Array elements contained solely in the x-z plane
%% Define Parameters
%Constants
c=11.81;%Gin/s
deg2rad=pi/180;
%Paramter Def’ns
alpha=90;%deg
alpharad=alpha*deg2rad;%rad
f=10;%GHz
lambda=c/f;%inches
k=2*pi/lambda;
theta=[-90:.1:90];
thetarad=theta*deg2rad;
thetao=0;
u=sin(thetarad);w=cos(thetarad);
uo=sin(thetao*deg2rad);
wo=cos(thetao*deg2rad);
rvec=[u;w];
R=11.46*lambda;

N=36;%Number of elements
d=1.5*lambda;
alphai=-0.5*alpharad + ([1:N] -1)*alpharad/(N-1);
xi=R*sin(alphai);
zi=R*cos(alphai);
denom=sqrt(xi.^2 + zi.^2);
nveci=[xi./denom; zi./denom];
EF=1.5;
%% Compute Conformal Pattern
elempat=(nveci.’*rvec).^(0.5*EF);
cosang=acos((nveci.’*rvec))/deg2rad;
indx=find(cosang > 90);
elempat(indx)=0;
phase=k*(xi’*(u-uo) + zi’*(w-wo));
Pat=sum(elempat.*exp(1i*phase));
[Pat_mat Pat_dB Pat_dBnorm]=process_vector(Pat);
%% Plot
figure 1),clf
set(gcf,’DefaultLineLineWidth’,2.5)
plot(theta,Pat_dBnorm,’-’,’color’,[0 0 1]),grid
axis([-90 90 -50 0])
set(gca,’FontSize’,16,’FontWeight’,’bold’)
set(gca,’XTick’,[-90:15:90])
title([‘Conformal Array Pattern’])
xlabel(‘\theta (degrees)’),ylabel(‘dB’)

