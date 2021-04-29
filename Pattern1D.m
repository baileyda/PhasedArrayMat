% 1D Pattern Code
% Computes Element Pattern (EP), Array Factor(AF)and array pattern (EP*AF)
% Arik D. Brown
 
clear all
 
%% Input Parameters
%ESA Parameters
 
%ESA opearating at tune freq
array_params.f=10;%Operating Frequency in GHz
array_params.fo=10;%Tune Frequency in GHz of the Phase Shifter,
 
array_params.nelem=30;%Number of Elements
array_params.d=0.5*(11.803/array_params.fo);%Element Spacing in Inches
 
array_params.EF=1.35;%EF
 
array_params.wgtflag=1;%0 = Uniform, 1 = Taylor Weighting
 
%$$$$These Parameters Only Used if array_params.wgtflag=1;
array_params.taylor.nbar=5;
array_params.taylor.SLL=30;%dB value
 
%Theta Angle Parameters
theta_angle.numpts=721;%Number of angle pts
theta_angle.min=-90;%degrees
theta_angle.max=90;%degrees
theta_angle.scan=0;%degrees
 
plotcommand.EP=0;%Plot EP if = 1
plotcommand.AF=0;%Plot EP if = 1
plotcommand.PAT=1;%Plot PAT if = 1
plotcommand.ALL=0;%Plot All patterns overlaid if = 1
 
%% Compute Patterns
 
if array_params.wgtflag==0
    array_params.amp_wgts=ones(array_params.nelem,1);
else
    array_params.amp_wgts=Taylor(array_params.nelem,array_params.taylor.SLL,...
                             array_params.taylor.nbar);    
end
 
theta_angle.vec=linspace(theta_angle.min,theta_angle.max,...
    theta_angle.numpts);%degrees
theta_angle.uvec=sind(theta_angle.vec);
theta_angle.uo=sind(theta_angle.scan);
 
%Initialize Element Pattern, Array Factor and Pattern
array.size=size(theta_angle.vec);
array.EP=zeros(array.size);%EP
array.AF=zeros(array.size);%AF
array.PAT=zeros(array.size);
 
%% Compute Patterns
 
%Compute AF1
[array.AF, array.AF_mag, array.AF_dB, array.AF_dBnorm]=...
    Compute_1D_AF(array_params.amp_wgts,array_params.nelem,...
    array_params.d,array_params.f,array_params.fo,...
    theta_angle.uvec,theta_angle.uo);
 
%Compute EP
[array.EP, array.EP_mag, array.EP_dB, array.EP_dBnorm]=...
    Compute_1D_EP(theta_angle.vec,array_params.EF);
 
%Compute PAT
[array.PAT, array.PAT_mag, array.PAT_dB, array.PAT_dBnorm] =...
    Compute_1D_PAT(array.EP,array.AF);
 
%% Plotting
 
if plotcommand.EP == 1
    %Plot EP in dB, Normalized
    figure,clf
    set(gcf,'DefaultLineLineWidth',2.5)
    plot(theta_angle.vec,array.EP_dBnorm,'--','color',[0 0 0]),hold
    grid
    axis([-90 90 -50 0])
    set(gca,'FontSize',16,'FontWeight','bold')
    title(['Element Pattern'])
    xlabel('\theta (degrees)'),ylabel('dB')
end
 
if plotcommand.AF == 1
    %Plot PAT in dB, Normalized
    figure,clf
    set(gcf,'DefaultLineLineWidth',2.5)
    plot(theta_angle.vec,array.AF_dBnorm,'color',[0 .7 0])
    grid
    axis([-90 90 -50 0])
    set(gca,'FontSize',16,'FontWeight','bold')
    title(['Linear ',num2str(array_params.nelem),' Element Array Array Factor'])
    xlabel('\theta (degrees)'),ylabel('dB')
end
 
if plotcommand.PAT == 1
    %Plot PAT in dB, Normalized
    figure,clf
    set(gcf,'DefaultLineLineWidth',2.5)
    plot(theta_angle.vec,array.PAT_dBnorm+array.EP_dBnorm,'color',[0 0 1]),hold
    grid
    axis([-90 90 -50 0])
    set(gca,'FontSize',16,'FontWeight','bold')
    title(['Linear ',num2str(array_params.nelem),' Element Array Pattern'])
    xlabel('\theta (degrees)'),ylabel('dB')
end
 
if plotcommand.ALL == 1
    %Plot ALL in dB, Normalized
    figure,clf
    set(gcf,'DefaultLineLineWidth',2.5)
    plot(theta_angle.vec,array.EP_dBnorm,'--','color',[0 0 0]),hold
    plot(theta_angle.vec,array.AF_dBnorm,'color',[0 .7 0])
    plot(theta_angle.vec,array.PAT_dBnorm+array.EP_dBnorm,'b-')
    grid
    axis([-90 90 -50 0])
%     axis([50 70 -20 0])
    set(gca,'FontSize',16,'FontWeight','bold')
    title(['Linear ',num2str(array_params.nelem),' Element Array'])
    xlabel('\theta (degrees)'),ylabel('dB')
    legend('EP','AF','PAT = EP * AF')
end

