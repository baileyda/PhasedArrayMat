% 1D Pattern Code
% Computes Patterns for Different Element Spacing to Illustrate Grating
% Lobes
% Arik D. Brown
clear all
%% Input Parameters
%ESA Parameters
%ESA opearating at tune freq
array_params.f=10;%Operating Frequency in GHz
array_params.fo=10;%Tune Frequency in GHz of the Phase Shifter,
array_params.nelem=30;%Number of Elements
array_params.d1=0.5*(11.803/array_params.fo);%Element Spacing in Inches
array_params.d2=1*(11.803/array_params.fo);%Element Spacing in Inches
array_params.d3=2*(11.803/array_params.fo);%Element Spacing in Inches
array_params.EF=1.35;%EF
array_params.amp_wgts=ones(array_params.nelem,1);
%Theta Angle Parameters
theta_angle.numpts=721;%Number of angle pts
theta_angle.min=-90;%degrees
theta_angle.max=90;%degrees
theta_angle.scan=0;%degrees
%% Compute Patterns
theta_angle.vec=linspace(theta_angle.min,theta_angle.max,...
theta_angle.numpts);%degrees
theta_angle.uvec=sind(theta_angle.vec);
theta_angle.uo=sind(theta_angle.scan);
%Initialize Element Pattern, Array Factor and Pattern
array.size=size(theta_angle.vec);
array.EP=zeros(array.size);%EP
array.AF1=zeros(array.size);%AF1
array.AF2=zeros(array.size);%AF2
array.AF3=zeros(array.size);%AF3
array.PAT1=zeros(array.size);%Pattern 1
array.PAT2=zeros(array.size);%Pattern 2
array.PAT3=zeros(array.size);%Pattern 3
%% Compute Patterns
%Compute AF1
[array.AF1, array.AF1_mag, array.AF1_dB, array.AF1_dBnorm]=...
Compute_1D_AF(array_params.amp_wgts,array_params.nelem,...
array_params.d1,array_params.f,array_params.fo,...
theta_angle.uvec,theta_angle.uo);
%Compute AF2
[array.AF2, array.AF2_mag, array.AF2_dB, array.AF2_dBnorm]=...
Compute_1D_AF(array_params.amp_wgts,array_params.nelem,...
array_params.d2,array_params.f,array_params.fo,...
theta_angle.uvec,theta_angle.uo);

%Compute AF3
[array.AF3, array.AF3_mag, array.AF3_dB, array.AF3_dBnorm]=...
Compute_1D_AF(array_params.amp_wgts,array_params.nelem,...
array_params.d3,array_params.f,array_params.fo,...
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
%Compute PAT3
[array.PAT3, array.PAT3_mag, array.PAT3_dB, array.PAT3_dBnorm] =...
Compute_1D_PAT(array.EP,array.AF3);
%% Plotting
%Plot PAT in dB, Normalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
plot(theta_angle.vec,array.PAT1_dBnorm+array.EP_dBnorm,’color’,[0 0 1]),hold
plot(theta_angle.vec,array.PAT2_dBnorm+array.EP_dBnorm,’color’,[0 .7 0]),
plot(theta_angle.vec,array.PAT3_dBnorm+array.EP_dBnorm,’color’,[1 0 0])
grid
axis([-90 90 -50 0])
set(gca,’FontSize’,16,’FontWeight’,’bold’)
title([‘Linear ‘,num2str(array_params.nelem),’ Element Array Pattern’])
xlabel(‘\theta (degrees)’),ylabel(‘dB’)
legend(‘d = 0.5*\lambda’,’d = 1*\lambda’,’d = 2*\lambda’)
