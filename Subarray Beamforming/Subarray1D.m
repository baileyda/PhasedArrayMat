%% Code to compute subarrayed architecture pattern
% Arik D. Brown
%% Input Parameters
IBWratio=1.1;%IBWratio -> f=fo
fo=4;%GHz Tune Frequency
f=IBWratio*fo;%GHz Operating Frequency
lambda=11.803/f;%inches
lambdao=11.803/fo;%inches
d=lambdao/2;%inches
theta=linspace(-90,90,721);%deg
thetao=30;%deg
u=sind(theta);
uo=sind(thetao);
SA.nelems=4;%Number of elements in Subarray
AF.nelems=4;%Number of elements in Backend AF
EF=1.5;
SA.wgts=ones(1,SA.nelems);
AF.wgts=ones(1,AF.nelems);
plotfigs.flags.EP=1;
plotfigs.flags.SA=1;
plotfigs.flags.AF=1;
plotfigs.flags.PAT=1;
plotfigs.flags.ALL=1;
plotfigs.axis.xlims=[-90 90];
plotfigs.axis.ylims=[-40 0];
%% Compute Pattern
% Element Pattern
[EP, EP_mag, EP_dB, EP_dBnorm] = Compute_1D_EP(theta,EF);
% Subarray AF
%!!! To simulate an array of elements without phase shifters input 0 for
%the last input to Compute_1D_AF
[SA.AF, SA.AF_mag, SA.AF_dB, SA.AF_dBnorm] =...
Compute_1D_AF(SA.wgts,SA.nelems,d,f,fo,u,uo);
%Backend AF
[AF.AF, AF.AF_mag, AF.AF_dB, AF.AF_dBnorm] =...
Compute_1D_AF(AF.wgts,AF.nelems,SA.nelems*d,f,f,u,uo);
%Pattern = Element Pattern x Subarray AF Pattern x AF Pattern
PAT=EP.*SA.AF.*AF.AF;
[PAT_mag,PAT_dB,PAT_dBnorm] = process_vector(PAT);
SA.scanvalue.afsa=SA.AF_dBnorm(u==uo);
EPscanvalue=EP_dBnorm(u==uo);
patnorm.sa=SA.scanvalue.afsa+EPscanvalue;
%% Plot Patterns
if plotfigs.flags.EP == 1
%Plot Pattern in dB, Unnormalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,EP_dBnorm),hold
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
t
itle([‘Element Pattern’],’FontSize’,14,’FontWeight’, ’bold’)
x
label(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.SA==1
%Plot Pattern in dB, Unnormalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,SA.AF_dBnorm),hold
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Subarray Pattern’],’FontSize’,14,’FontWeight’, ’bold’)
x
label(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.AF==1
%Plot Pattern in dB, Unnormalized
figure,clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,AF.AF_dBnorm),hold
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Array Pattern’],’FontSize’,14,’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.PAT==1
%Plot Pattern in dB, Unnormalized
figure 4),clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,PAT_dBnorm+patnorm.sa),hold
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
t
itle([‘Array Pattern’],’FontSize’,14,’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
set(gca,’Position’,[0.13 0.11 0.775 0.815])
end
if plotfigs.flags.ALL==1
%Plot Pattern in dB, Unnormalized
figure 5),clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’,’bold’)
plot(theta,EP_dBnorm,’color’,[0 0 0]),hold
plot(theta,SA.AF_dBnorm,’color’,[0 .7 0])
plot(theta,AF.AF_dBnorm,’color’,[.7 0 1])
plot(theta,PAT_dBnorm+patnorm.sa,’color’,[0 0 1])
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Composite Array Pattern, f = ‘,num2str(IBWratio),’*f_{o}’],...
‘FontSize’,14,’FontWeight’,’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
legend(‘EP’,’SA.AF’,’AF’,’Pattern=EP*SA.AF*AF’)
end
