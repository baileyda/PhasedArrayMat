%% Code to compute and plot grating lobes in sine space
%Arik D. Brown
%% Define input parameters
c=11.803;%Gin/s
f=3;%GHz
lambda=c/f;%in

%Element Spacing
dx=.5*lambda;%in
dy=.29*lambda;%in
%Grid Type
elemgrid=2;%1 = Rectangular Grid, 2 = Triangular Grid
%Define Scan Angle
thetaELo=0;%deg thetaEL
thetaAZo=0;%deg phiAZ
mainbeam.uo=cosd(thetaELo)*sind(thetaAZo);
mainbeam.vo=sind(thetaELo);
%Define Unit Circle
phi=[0:1:360];%degrees
unitcircle.u=cosd(phi);
unitcircle.v=sind(phi);
%Define Main Beam Unit Circle
mainbeam.u=mainbeam.uo + unitcircle.u;
mainbeam.v=mainbeam.vo + unitcircle.v;
%Define Grating Lobe Locations
if elemgrid == 1
gratinglobe.uo=mainbeam.uo+[ 0 -1 -1 1 1 0 -1 1]*lambda/dx;
gratinglobe.vo=mainbeam.vo+[ 1 -1 0 0 1 -1 1 -1]*lambda/dy;
plotfigs.titletxt=’Rectangular Grid’;
plotfigs.markertype=’s’;
elseif elemgrid ==2
gratinglobe.uo=mainbeam.uo+[ 0 0 1 -1 -1 1 2 -2]*lambda/(2*dx);
gratinglobe.vo=mainbeam.vo+[ 2 -2 1 -1 1 -1 0 0]*lambda/(2*dy);
plotfigs.titletxt=’Triangular Grid’;
plotfigs.markertype=’^’;
end
%% Plot Grating Lobes
plotfigs.limit=4.5;
figure 1)
clf
%Mainbeam Unit Circle
plot(mainbeam.u,mainbeam.v,’-’,’color’,[0 0 0],’LineWidth’,1.5),hold

%Main Beam
plot(mainbeam.uo,mainbeam.vo,plotfigs.markertype,’color’,[0 0 0],’MarkerSize’,6,’LineWidth’,1.5)
%Grating Lobes
plot(gratinglobe.uo,gratinglobe.vo,plotfigs.markertype,’color’,[1 0 0],’LineWidth’,1.5)
%Grating Lobe Unit Circles
for ig=1:length(gratinglobe.uo)
if elemgrid == 1
plot(unitcircle.u+ gratinglobe.uo(ig),unitcircle.v+gratinglobe.vo(ig),’r-’,’LineWidth’,1.5)
elseif elemgrid == 2
plot(unitcircle.u+ gratinglobe.uo(ig),unitcircle.v+gratinglobe.vo(ig),’r-’,’LineWidth’,1.5)
end
end
axis([-plotfigs.limit plotfigs.limit -plotfigs.limit plotfigs.limit])
set(gca,’XTick’,[-plotfigs.limit:.5:plotfigs.limit],’YTick’,[-plotfigs.limit:.5:plotfigs.limit])
xlabel(‘u’,’FontSize’,14,’FontWeight’,’bold’)
ylabel(‘v’,’FontSize’,14,’FontWeight’,’bold’)
title(plotfigs.titletxt,’FontSize’,18,’FontWeight’,’bold’)
grid
set(gca,’Fontsize’,14,’Fontweight’,’bold’,’linewidth’,1.0)

