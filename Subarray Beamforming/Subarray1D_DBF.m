%% Code to compute subarrayed architecture pattern
% Arik D. Brown
%% Input Parameters
IBWratio=1.;%IBWratio -> f=fo
fo=4;%GHz Tune Frequency
f=IBWratio*fo;%GHz Operating Frequency
lambda=11.803/f;%inches
lambdao=11.803/fo;%inches
d=lambdao/2;%inches
theta=linspace(-90,90,721);%deg
thetao1=-6;%deg
thetao2=0;%deg
thetao3=6;%deg
u=sind(theta);
uo1=sind(thetao1);
uo2=sind(thetao2);
uo3=sind(thetao3);
SA.nelems=4;%Number of elements in Subarray
AF.nelems=4;%Number of elements in Backend AF
EF=1.5;
SA.wgts=ones(1,SA.nelems);
AF.wgts=ones(1,AF.nelems);
plotfigs.flags.EP=1;
plotfigs.flags.SA=1;
plotfigs.flags.AF=1;
plotfigs.flags.PATs=1;
plotfigs.axis.xlims=[-90 90];
plotfigs.axis.ylims=[-40 0];
%% Compute Pattern
% Element Pattern
[EP, EP_mag, EP_dB, EP_dBnorm] = Compute_1D_EP(theta,EF);
% Subarray AF
%!!! To simulate an array of elements without phase shifters input 0 for
%the last input to Compute_1D_AF
[SA.AF, SA.AF_mag, SA.AF_dB, SA.AF_dBnorm] =...
Compute_1D_AF(SA.wgts,SA.nelems,d,f,fo,u,uo2);
%Backend AFs for DBF
%%Beam 1
[AF.AF1, AF.AF1_mag, AF.AF1_dB, AF.AF1_dBnorm] =...
Compute_1D_AF(AF.wgts,AF.nelems,SA.nelems*d,f,f,u,uo1);
%%Beam 2
[AF.AF2, AF.AF2_mag, AF.AF2_dB, AF.AF2_dBnorm] =...
Compute_1D_AF(AF.wgts,AF.nelems,SA.nelems*d,f,f,u,uo2);
%%Beam 3
[AF.AF3, AF.AF3_mag, AF.AF3_dB, AF.AF3_dBnorm] =...
Compute_1D_AF(AF.wgts,AF.nelems,SA.nelems*d,f,f,u,uo3);
%Pattern = Element Pattern x Subarray AF Pattern x AF Pattern
PAT1=EP.*SA.AF.*AF.AF1;
[PAT1_mag,PAT1_dB,PAT1_dBnorm] = process_vector(PAT1);
PAT2=EP.*SA.AF.*AF.AF2;
[PAT2_mag,PAT2_dB,PAT2_dBnorm] = process_vector(PAT2);
PAT3=EP.*SA.AF.*AF.AF3;
[PAT3_mag,PAT3_dB,PAT3_dBnorm] = process_vector(PAT3);
SA.scanvalue.afsa1=SA.AF_dBnorm(u==uo1);
SA.scanvalue.afsa2=SA.AF_dBnorm(u==uo2);
SA.scanvalue.afsa3=SA.AF_dBnorm(u==uo3);
EPscanvalue=EP_dBnorm(u==uo);
patnorm.sa1=SA.scanvalue.afsa1+EPscanvalue;
patnorm.sa2=SA.scanvalue.afsa2+EPscanvalue;
patnorm.sa3=SA.scanvalue.afsa3+EPscanvalue;
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
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’,’bold’)
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
t
itle([‘Array Pattern’],’FontSize’,14,’FontWeight’, ’bold’)
xlabel(‘\theta (degrees)’,’FontSize’,12,’FontWeight’, ’bold’)
ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
end
if plotfigs.flags.PATs==1
%Plot Pattern in dB, Unnormalized
figure 5),clf
set(gcf,’DefaultLineLineWidth’,1.5)
s
et(gcf,’DefaultTextFontSize’,12,’DefaultTextFontWeight’, ’bold’)
plot(theta,SA.AF_dBnorm,’color’,[0 0 0]),hold
plot(theta,PAT1_dBnorm+patnorm.sa1,’--’,’color’,[.7 0 1])
plot(theta,PAT2_dBnorm+patnorm.sa2,’-’,’color’,[0 0 1])
plot(theta,PAT3_dBnorm+patnorm.sa3,’--’,’color’,[0 .7 0])
grid
axis([plotfigs.axis.xlims plotfigs.axis.ylims])
title([‘Composite Array Pattern’],’FontSize’,14, ’FontWeight’, ’bold’)
x
label(‘\theta (degrees)’,’FontSize’,12,’FontWeight’, ’bold’)

ylabel(‘dB’,’FontSize’,12,’FontWeight’,’bold’)
set(gca,’FontSize’,12,’FontWeight’,’bold’)
set(gcf,’color’,’white’)
legend(‘SA.AF’,’Beam1’,’Beam2’,’Beam3’)
end