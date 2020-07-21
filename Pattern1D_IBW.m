% 1D Pattern Code Demonstrating Beamsquint Due to IBW Constraints
% Arik D. Brown

clear all
%% Input Parameters
%ESA Parameters
%ESA opearating at tune freq
array_params.f1=10;%Operating Frequency in GHz
array_params.deltaf=-0.2;%
array_params.f2=10+array_params.deltaf;%Operating Frequency in GHz (Squinted Beam)
array_params.fo=10;%Tune Frequency in GHz of the Phase Shifter,
array_params.nelem=30;%Number of Elements
array_params.d=0.5*(11.803/array_params.fo);%Element Spacing in Inches
array_params.EF=1.35;%EF
array_params.select_wgts=0;
array_params.amp_wgts=ones(array_params.nelem,1);
%Theta Angle Parameters
theta_angle.numpts=1001;%Number of angle pts
theta_angle.min=10;%degrees
theta_angle.max=50;%degrees
theta_angle.scan=30;%degrees
%% Compute Patterns
theta_angle.vec=linspace(theta_angle.min,theta_angle.max,...
theta_angle.numpts);%degrees
theta_angle.uvec=sind(theta_angle.vec);
theta_angle.uo=sind(theta_angle.scan);
%Initialize Element Pattern, Array Factor and Pattern
array.size=size(theta_angle.vec);
array.EP=zeros(array.size);%EP
array.AF1=zeros(array.size);%AF1 f=fo
array.AF2=zeros(array.size);%AF2 f=fo+deltaf
array.PAT=zeros(array.size);
%% Compute Patterns
%Compute AF1
[array.AF1, array.AF1_mag, array.AF1_dB, array.AF1_dBnorm]=...
Compute_1D_AF(array_params.amp_wgts,array_params.nelem,...
array_params.d,array_params.f1,array_params.fo,...
theta_angle.uvec,theta_angle.uo);

%Compute AF2
[array.AF2, array.AF2_mag, array.AF2_dB, array.AF2_dBnorm]=...
Compute_1D_AF(array_params.amp_wgts,array_params.nelem,...
array_params.d,array_params.f2,array_params.fo,...
theta_angle.uvec,theta_angle.uo);
%Compute EP
[array.EP, array.EP_mag, array.EP_dB, array.EP_dBnorm]=...
Compute_1D_EP(theta_angle.vec,array_params.EF);
%Compute PAT1
[array.PAT1, array.PAT1_mag, array.PAT1_dB, array.PAT1_dBnorm] =...
Compute_1D_PAT(array.EP,array.AF1);
%Compute PAT2
[array.PAT2, array.PAT2_mag, array.PAT2_dB, array.PAT2_dBnorm] =...
Compute_1D_PAT(array.EP,array.AF2);
%% Plotting
%Plot PAT1 and PAT2 in dB, Normalized
figure 1),clf
set(gcf,’DefaultLineLineWidth’,2.5)
plot(theta_angle.vec,array.PAT1_dBnorm,’k-’),hold
plot(theta_angle.vec,array.PAT2_dB - max(array.PAT1_dB),’k--’),hold
grid
axis([25 35 -5 0])
set(gca,’FontSize’,16,’FontWeight’,’bold’)
title([‘Linear ‘,num2str(array_params.nelem),’ Element Array’])
xlabel(‘\theta (degrees)’),ylabel(‘dB’)
legend(‘f = f_{o}’,’f = f_{o} + \Delta f’)
set(gca,’XTick’,[25:1:35],’YTick’,[-5:1:0])
