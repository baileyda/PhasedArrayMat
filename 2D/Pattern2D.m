% 2D Pattern Code
% Computes Element Pattern (EP), Array Factor(AF)and array pattern (EP*AF)
% Arik D. Brown
clear all
%% Input Parameters
%ESA Parameters
%ESA opearating at tune freq
array_params.f=3;%Operating Frequency in GHz
array_params.fo=3;%Tune Frequency in GHz of the Phase Shifter,
array_params.nelem.x=20;%Number of Elements in x
array_params.nelem.y=20;%Number of Elements in y
array_params.d.x=(1/(1+sind(90)))*(11.803/array_params.fo);%Element Spacing in Inches
array_params.d.y=(1/(1+sind(90)))*(11.803/array_params.fo);%Element Spacing in Inches
array_params.EF=1.5;%EF
array_params.flag.gain=0;%0 = Don’t Compute Gain, 1=Compute Gain
array_params.flag.wgt=0;%0 = Uniform, 1 = Taylor Weighting
array_params.flag.error_rand=0;%0 = Ideal (no errors), 1 = Random Phase/Amp. errors
array_params.flag.error_quant=0;%0 = Ideal (no quant.), 1 = Quant.
array_params.tilt.option=0;%0 = Ideal position (No roll, pitch or yaw)
%1 = Roll (rotation about z axis)
%2 = Pitch (rotation about x axis)
%3 = Yaw (rotation about y axis)
array_params.tilt.angle=30;%degrees
array_params.error.rand.amp.onesigma=0.5;%dB
array_params.error.rand.phs.onesigma=6;%degrees
array_params.bits=6;
%$$$$These Parameters Only Used if array_params.wgtflag=1;
array_params.taylor.nbar=5;
array_params.taylor.SLL=30;%dB value
%Theta and Phi Angle Parameters (Antenna Coordinates)
theta_angle.numpts=361;%Number of angle pts in theta
phi_angle.numpts=361;%Number of angle pts in phi
theta_angle.min=0;%degrees
theta_angle.max=90;%degrees
phi_angle.min=0;%degrees
phi_angle.max=360;%degrees
theta_angle.scan=45;%degrees
phi_angle.scan=90;%degrees
plotcommand.coord=0;%0 = Antenna Coord., 1 = Radar Coordinates
plotcommand.error=0;%0 = Don’t Plot Errors, 1 = Plot Errors
plotcommand.EP=0;%Plot EP if = 1
plotcommand.AF=0;%Plot AF if = 1
plotcommand.PAT=1;%Plot PAT if = 1
plotcommand.INTGAIN=0;%Plot INTGAIN if = 1
array_params.dBfloor.EP=-20;% dB value for plotting
array_params.dBfloor.PAT=-50;% dB value for plotting

%% Computations

if array_params.flag.wgt==0
array_params.amp_wgts.mat=ones(array_params.nelem.y,array_params.nelem.x);
else
array_params.amp_wgts.vec.x=Taylor(array_params.nelem.x,array_params.taylor.SLL,...
array_params.taylor.nbar);
array_params.amp_wgts.vec.y=Taylor(array_params.nelem.y,array_params.taylor.SLL,...
array_params.taylor.nbar);
array_params.amp_wgts.mat=array_params.amp_wgts.vec.y*...
array_params.amp_wgts.vec.x’;
end
if array_params.flag.error_rand==0
array_params.error.rand.amp.mat=ones(size(array_params.amp_wgts.mat));

array_params.error.rand.amp.mat_dB=10*log10(array_params.error.rand.amp.mat);%dB
array_params.error.rand.phs.mat=zeros(size(array_params.amp_wgts.mat));%degrees
elseif array_params.flag.error_rand==1
array_params.error.rand.amp.mat_dB=...
array_params.error.rand.amp.onesigma*randn(size(array_params.amp_wgts.mat));%dB
array_params.error.rand.amp.mat=...
10.^(array_params.error.rand.amp.mat_dB/10);
array_params.error.rand.phs.mat=...
array_params.error.rand.phs.onesigma*randn(size(array_params.amp_wgts.mat));%degrees
end
theta_angle.vec=linspace(theta_angle.min,theta_angle.max,...
theta_angle.numpts);%degrees
phi_angle.vec=linspace(phi_angle.min,phi_angle.max,...
phi_angle.numpts);%degrees
[theta_angle.mat phi_angle.mat]=meshgrid(theta_angle.vec,phi_angle.vec);
if array_params.tilt.option == 0
array_params.tilt_mat=[1 0 0;...
0 1 0;...
0 0 1];
elseif array_params.tilt.option == 1
array_params.tilt_mat=[cosd(array_params.tilt.angle) -sind(array_params.tilt.angle) 0;...
sind(array_params.tilt.angle) cosd(array_params.tilt.angle) 0;...
0 0 1];
elseif array_params.tilt.option == 2
array_params.tilt_mat=[1 0 0;...
0 cosd(array_params.tilt.angle) sind(array_params.tilt.angle);...
0 -sind(array_params.tilt.angle) cosd(array_params.tilt.
angle)];
elseif array_params.tilt.option == 3
array_params.tilt_mat=[cosd(array_params.tilt.angle) 0 -sind(array_params.tilt.angle);...
0 1 0;...

sind(array_params.tilt.angle) 0 cosd(array_params.tilt.angle)];
end
sinespace.umat=sind(theta_angle.mat).*cosd(phi_angle.mat);
sinespace.vmat=sind(theta_angle.mat).*sind(phi_angle.mat);
sinespace.wmat=sqrt(abs(1-(sinespace.umat.^2)-(sinespace.vmat.^2)));
sinespace.uo=sind(theta_angle.scan)*cosd(phi_angle.scan);
sinespace.vo=sind(theta_angle.scan)*sind(phi_angle.scan);
sinespace.wo=sqrt(1-(sinespace.uo.^2)-(sinespace.vo.^2));
sinespace.uvwmat=[sinespace.umat(:).’;sinespace.vmat(:).’;...
sinespace.wmat(:).’];
sinespace.uvwmatnew=array_params.tilt_mat*sinespace.uvwmat;
sinespace.umatnew=reshape(sinespace.uvwmatnew(1,:),size(sinespace.umat));
sinespace.vmatnew=reshape(sinespace.uvwmatnew(2,:),size(sinespace.vmat));
sinespace.wmatnew=reshape(sinespace.uvwmatnew(3,:),size(sinespace.wmat));
if plotcommand.coord==0
plotxy.x=sinespace.umatnew;
plotxy.y=sinespace.vmatnew;
plotxy.xtxt=’u’;
plotxy.ytxt=’v’;
plotxy.xtick.min=-1;
plotxy.xtick.max=1;
plotxy.ytick.min=-1;
plotxy.ytick.max=1;
plotxy.tick.delta=.25;
elseif plotcommand.coord==1
radarcoord.thetaAZmat=atan2(sinespace.umat,sinespace.wmat)*180/pi;%degrees
radarcoord.thetaELmat=asind(sinespace.vmat);%degrees
plotxy.x=radarcoord.thetaAZmat;
plotxy.y=radarcoord.thetaELmat;
plotxy.xtxt=’\theta_{AZ}’;
plotxy.ytxt=’\theta_{EL}’;
plotxy.xtick.min=-90;
plotxy.xtick.max=90;
plotxy.ytick.min=-90;
plotxy.ytick.max=90;
plotxy.tick.delta=30;
end

%Initialize Element Pattern, Array Factor and Pattern
array.size=size(theta_angle.mat);
array.EP=zeros(array.size);%EP
array.AF=zeros(array.size);%AF
array.PAT=zeros(array.size);
%% Compute Patterns
%Compute AF1
if array_params.flag.error_quant == 0
[array.AF, array.AF_mag, array.AF_dB, array.AF_dBnorm]=...
Compute_2D_AF(array_params.amp_wgts.mat,...
array_params.nelem.x,array_params.nelem.y,...
array_params.d.x,array_params.d.y,...
array_params.f,array_params.fo,...
array_params.error.rand.amp.mat,array_params.error.rand.phs.mat,...
sinespace.umat,sinespace.vmat,...
sinespace.uo,sinespace.vo);
elseif array_params.flag.error_quant == 1
[array.AF, array.AF_mag, array.AF_dB, array.AF_dBnorm]=...
Compute_2D_AF_quant(array_params.amp_wgts.mat,...
array_params.nelem.x,array_params.nelem.y,...
array_params.d.x,array_params.d.y,...
array_params.f,array_params.fo,...
array_params.error.rand.amp.mat,array_params.error.rand.phs.mat,...
array_params.bits,...
sinespace.umat,sinespace.vmat,...
sinespace.uo,sinespace.vo);
end
%Compute EP
[array.EP, array.EP_mag, array.EP_dB, array.EP_dBnorm]=...
Compute_2D_EP(theta_angle.mat,array_params.EF);
%Compute PAT
[array.PAT, array.PAT_mag, array.PAT_dB, array.PAT_dBnorm] =...
Compute_2D_PAT(array.EP,array.AF);
if array_params.flag.gain==1
[array.INTGAINpeak array.IdealGain array.INTGAINpeakdB array.IdealGaindB...
array.INTGAIN array.INTGAIN_dB array.INTGAIN_dBnorm] =...
Compute_2D_INTGAIN(array.PAT_mag,...
theta_angle.vec,theta_angle.vec,theta_angle.mat,...
array_params.nelem.x,array_params.nelem.y,...

array_params.d.x,array_params.d.y,...
array_params.f);
[array.INTGAINpeakdB array.IdealGaindB]
end
%% Plotting
% close all
if plotcommand.error == 1
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
imagesc(array_params.error.rand.amp.mat_dB)
shading interp
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Random Amplitude Error(dB)’)
xlabel(‘Elements in x’),ylabel(‘Elements in y’)
view(2)
colorbar
caxis([-2 2])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
set(gcf, ‘color’, ‘white’);
h=figure;clf
binvector=[-2:.1:2];
set(gcf,’DefaultLineLineWidth’,1.5)
hist(array_params.error.rand.amp.mat_dB(:),binvector);
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Amplitude Error Distribution’)
xlabel(‘Amplitude Error Bins (dB)’)
axis tight
set(gca,’XTick’,[-2:0.5:2])
set(gcf, ‘color’, ‘white’);
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
imagesc(array_params.error.rand.phs.mat)
shading interp
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Random Phase Error(^{o})’)
xlabel(‘Elements in x’),ylabel(‘Elements in y’)
view(2)
colorbar
caxis([-20 20])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
set(gcf, ‘color’, ‘white’);

h=figure;clf
binvector=[-20:1:20];
set(gcf,’DefaultLineLineWidth’,1.5)
hist(array_params.error.rand.phs.mat(:),binvector);
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Phase Error Distribution’)
xlabel(‘Phase Error Bins (^{o})’)
axis tight
set(gca,’XTick’,[-20:5:20])
set(gcf, ‘color’, ‘white’);
end
if plotcommand.EP == 1
%Plot EP in dB, Normalized
plotEP=array.EP_dBnorm;
plotEP(array.EP_dBnorm < array_params.dBfloor.PAT)=array_params.dBfloor.PAT;
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,array.EP_dBnorm),hold
shading interp
colorbar
caxis([array_params.dBfloor.EP 0])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Element Pattern’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(2)%Plot EP in dB, Normalized
axis tight
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
set(gcf, ‘color’, ‘white’);
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,plotEP),hold
shading interp
colorbar
caxis([array_params.dBfloor.PAT 0])
zlim([array_params.dBfloor.PAT 0])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Element Pattern’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(3)

set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
set(gcf, ‘color’, ‘white’);
end
if plotcommand.AF == 1
%Plot PAT in dB, Normalized
plotAF=array.AF_dBnorm;
plotAF(array.AF_dBnorm < array_params.dBfloor.PAT)=array_params.dBfloor.PAT;
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,array.AF_dBnorm)
shading interp
colorbar
caxis([array_params.dBfloor.PAT 0])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘AF’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(2)
axis tight
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
set(gcf, ‘color’, ‘white’);
%Plot PAT in dB, Normalized
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,plotAF)
shading interp
colorbar
caxis([array_params.dBfloor.PAT 0])
zlim([array_params.dBfloor.PAT 0])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘AF’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(3)
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
set(gcf, ‘color’, ‘white’);
end

if plotcommand.PAT == 1
%Plot PAT in dB, Normalized
plotPAT=array.PAT_dBnorm;
plotPAT(array.PAT_dBnorm <= array_params.dBfloor.PAT)=array_params.dBfloor.PAT;
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,array.PAT_dBnorm)
shading interp
colorbar
caxis([array_params.dBfloor.PAT 0])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Pattern’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(2)
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.m axis tight
set(gcf, ‘color’, ‘white’);
%Plot PAT in dB, Normalized
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,plotPAT)
shading interp
colorbar
caxis([array_params.dBfloor.PAT 0])
zlim([array_params.dBfloor.PAT 0])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Pattern’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(3)
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
set(gcf, ‘color’, ‘white’);
end
if array_params.flag.gain == 1 && plotcommand.INTGAIN == 1
%Plot INTGAIN in dB, Normalized
plotINTGAIN=array.INTGAIN_dB;
plotINTGAIN(array.INTGAIN_dB <= array.INTGAINpeakdB+array_params.dBfloor.PAT)=...
array.INTGAINpeakdB+array_params.dBfloor.PAT;

h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,array.INTGAIN_dB)
shading interp
colorbar
caxis([array.INTGAINpeakdB+array_params.dBfloor.PAT array.INTGAINpeakdB])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Integrated Gain’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(2)
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
axis tight
set(gcf, ‘color’, ‘white’);
%Plot INTGAIN in dB, Normalized
h=figure;clf
set(gcf,’DefaultLineLineWidth’,1.5)
surf(plotxy.x,plotxy.y,plotINTGAIN)
shading interp
colorbar
caxis([array.INTGAINpeakdB+array_params.dBfloor.PAT array.INTGAINpeakdB])
zlim([array.INTGAINpeakdB+array_params.dBfloor.PAT array.INTGAINpeakdB])
set(gca,’FontSize’,14,’FontWeight’,’bold’)
title(‘Integrated Gain’)
xlabel(plotxy.xtxt),ylabel(plotxy.ytxt),zlabel(‘dB’)
view(3)
set(gca,’XTick’,[plotxy.xtick.min:plotxy.tick.delta:plotxy.xtick.max])
set(gca,’YTick’,[plotxy.ytick.min:plotxy.tick.delta:plotxy.ytick.max])
set(gcf, ‘color’, ‘white’);
end




