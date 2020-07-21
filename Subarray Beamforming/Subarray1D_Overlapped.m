%% Code to compute subarrayed architecture pattern
% Arik D. Brown
%% Input Parameters
IBWratio=1.1;%IBWratio -> f=fo
fo=4;%GHz Tune Frequency
f=IBWratio*fo;%GHz Operating Frequency
lambda=11.803/f;%inches
lambdao=11.803/fo;%inches
d=lambdao/2;%inches
theta=linspace(-10,50,721);%deg
thetao=30;%deg
u=sind(theta);
uo=sind(thetao);
OSA.ratio=2;
SA.nelems=8;%Number of elements in Subarray
SA.nsas=20;%Number of Subarrays combined in Backend AF
O SA.nelems=OSA.ratio*SA.nelems;%Number of elements in Overlapped Subarray(OSA)
O SA.nsas=SA.nsas-(OSA.ratio-1);%Number of OSAs combined in Backend OSA AF
EF=1*1.5;
SA.wgts.elems=ones(1,SA.nelems);
% SA.wgts.elems=Taylor(SA.nelems,35,6);
SA.wgts.sas=ones(1,SA.nsas);
% OSA.wgts.elems=ones(1,OSA.nelems);
OSA.wgts.elems=Taylor(OSA.nelems,35,6);
OSA.wgts.sas=ones(1,OSA.nsas);
plotfigs.flags.EP=1;
plotfigs.flags.SA=1;
plotfigs.flags.AF=1;
plotfigs.flags.PATs=1;
plotfigs.flags.ALL=1;
plotfigs.axis.xlims=[min(theta) max(theta)];
plotfigs.axis.ylims=[-60 0];
%% Compute Pattern
% Element Pattern
[EP, EP_mag, EP_dB, EP_dBnorm] = Compute_1D_EP(theta,EF);
% Subarray and Overlapped Subarrays
% Subarray AF
%!!! To simulate an array of elements without phase shifters input 0 for
%the last input to Compute_1D_AF
[SA.AFsa, SA.AFsa_mag, SA.AFsa_dB, SA.AFsa_dBnorm] =...
Compute_1D_AF(SA.wgts.elems,SA.nelems,d,f,fo,u,uo);
%OSA AF
[OSA.AFsa, OSA.AFsa_mag, OSA.AFsa_dB, OSA.AFsa_dBnorm] =...
Compute_1D_AF(OSA.wgts.elems,OSA.nelems,d,f,fo,u,uo);
%Backend AFs for Subarray and OSA
%Subarray Beamforming
[SA.AF, SA.AF_mag, SA.AF_dB, SA.AF_dBnorm] =...
Compute_1D_AF(SA.wgts.sas,SA.nsas,SA.nelems*d,f,f,u,uo);
%Overlapped Subarray Beamforming
[OSA.AF, OSA.AF_mag, OSA.AF_dB, OSA.AF_dBnorm] =...
Compute_1D_AF(OSA.wgts.sas,OSA.nsas,SA.nelems*d,f,f,u,uo);
%Pattern = Element Pattern x Subarray AF Pattern x AF Pattern
SA.PAT=EP.*SA.AFsa.*SA.AF;
[SA.PAT_mag,SA.PAT_dB,SA.PAT_dBnorm] = process_vector(SA.PAT);
OSA.PAT=EP.*OSA.AFsa.*OSA.AF;
[ OSA.PAT_mag,OSA.PAT_dB,OSA.PAT_dBnorm] = process_vector(OSA.PAT);
SA.scanvalue.afsa=SA.AFsa_dBnorm(u==uo);
OSA.scanvalue.afsa=OSA.AFsa_dBnorm(u==uo);
EPscanvalue=EP_dBnorm(u==uo);
patnorm.sa=SA.scanvalue.afsa+EPscanvalue;
patnorm.osa=OSA.scanvalue.afsa+EPscanvalue;
%% Plot Patterns
if plotfigs.flags.EP == 1
%Plot Pattern in dB, Normalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,EP_dBnorm),hold
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Element Pattern’],’FontSize’,14,’FontWeight’, ’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’, ’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.SA==1
%Plot SA Pattern in dB, Normalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,SA.AFsa_dBnorm),hold
plot(theta,OSA.AFsa_dBnorm,’--’)
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Subarray and OSA AF Patterns’],’FontSize’,14, ’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’, ’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.AF==1
%Plot AF in dB, Normalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,SA.AF_dBnorm),hold
plot(theta,OSA.AF_dBnorm,’--’)
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘AF Pattern’],’FontSize’,14,’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,EP_dBnorm),hold
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Element Pattern’],’FontSize’,14,’FontWeight’, ’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’, ’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.SA==1
%Plot SA Pattern in dB, Normalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,SA.AFsa_dBnorm),hold
plot(theta,OSA.AFsa_dBnorm,’--’)
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Subarray and OSA AF Patterns’],’FontSize’,14, ’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’, ’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.AF==1
%Plot AF in dB, Normalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,SA.AF_dBnorm),hold
plot(theta,OSA.AF_dBnorm,’--’)
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘AF Pattern’],’FontSize’,14,’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end

set(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,SA.AFsa_dBnorm,’k-’),hold
plot(theta,SA.AF_dBnorm,’-’,’color’,[0 0.7 0])
plot(theta,SA.PAT_dBnorm+patnorm.sa,’-’,’color’,[0 0 1])
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Composite Array Pattern with No Overlap and f = ‘,num2str(IBWratio),’*f_{o}’],...
‘FontSize’,14,’FontWeight’,’bold’)
x
label(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
legend(‘SA Subarray Pattern’,’SA Backend AF Pattern’,’SA Total Pattern’)
end




