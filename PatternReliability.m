%Pattern code for an mXn Array
%Element and Subarray Failures can be evaluated
%Written by: Jabberia Miller
clear all;clc
%% Variables
f=3.85e9; % Frequency (GHz)
c=11.80285e9; %Speed of Light (in/s)
lambda=c/f;
k=2*pi/lambda; %wavenumber
omega=2*pi*f;
elemfac=1.35; %Element Factor
theta_scan = 0; %Theta Scan Angle (degrees)
phi_scan = 0; %Phi Scan Angle (degrees)
%Number of Elements in x and y
M=40; %number of x elements
N=40; %number of y elements
subM = 5; %number of subarrays in x - must be exactly divisible by num of elements
subN = 5; %number of subarrays in y - must be exactly divisible by num of elements
removeElem = 2; % 0 - no failures
% 1 - fail random elements
% 2 - fail random subarray
percentFail=10; % percent of elements to fail
%%
%Compute Spacing between elements
dx=lambda/2;
dy=lambda/2;
%Compute Max Gain
A=(M*dy)*(N*dx);
G=4*pi*A/(lambda^2);
%Grid of locations for element x and y
xtemp = ([1:M] - 0.5*(M+1))*dx;
ytemp = ([1:N] - 0.5*(N+1))*dy;
[x y] = meshgrid(xtemp,ytemp);
%Grid of locations for subarray x and y
subdx = M/subM ; %number of elements per subarray in x
subdy = N/subN ; %number of elements per subarray in y
subtempx = ([1:subM] - 0.5*(subM+1))*subdx*dx;
subtempy = ([1:subN] - 0.5*(subN+1))*subdy*dy;
[subx suby] = meshgrid(subtempx,subtempy);
%Associate Elements with subarrays
subnum = subM*subN; %Calculate number of Subarrays
sub = ones(size(reshape(x,M*N,1))); %storage variable
x1=reshape(x,M*N,1); % reshape x array into vector
y1=reshape(y,M*N,1); % reshape y array into vector
offst=0.1;
for ict = 1:subnum
jj = find(x1>=(subx(ict)-subdx*dx/2-offst) & x1<=(subx(ict)+subdx*dx/2+offst) & y1 <=(suby(ict)+subdy*dy/2+offst) & y1 >=(suby(ict)-subdy*dy/2-offst));
sub(jj) = ict;
end
sub = reshape(sub,M,N);
%compute angles where pattern will be evaluated over
[Az,El] = meshgrid(-90:0.2:90,-90:0.2:90);
%convert to Az/EL angles to theta/phi
Phyr = Az * pi/180; %convert to radians
Thyr = (90 - El) * pi/180; %convert to radians
u = sin(Thyr) .* cos(Phyr);
v = sin(Thyr) .* sin(Phyr);
w = cos(Thyr);
thetamat = acos(u ./ sqrt(v.^2 + w.^2 + u.^2)) ;
phimat = atan2(w,v) ;
%Compute direction cosines
u=sin(thetamat).*cos(phimat);
v=sin(thetamat).*sin(phimat);
%Compute steering angles
thetao=theta_scan*pi/180;
phio=-phi_scan*pi/180;
uo=sin(thetao)*cos(phio);
vo=sin(thetao)*sin(phio);
%Remove elements or subarrays
AF=zeros(size(u));
elempat=zeros(size(u));
wgts = ones(size(x));
switch removeElem
case 1 %remove random elements
numdelete = round(percentFail*M*N/100);
DeleteList = randperm(M*N);
DeleteList = DeleteList(1:numdelete);
wgts(DeleteList)= 0;
case 2 % remove subarray
deletesub = randperm(subM*subN);
[tmpx tmpy] = find(sub == deletesub(1));
wgts(tmpx,tmpy)=0;
end
%Compute Array Factor
for ii=1:M*N
tm = (x(ii)*u + y(ii)*v)/c;
tmo = (x(ii)*uo + y(ii)*vo)/c;
AF = AF + wgts(ii)*exp(j*omega*(tm-tmo));
end
%Compute element pattern
elempat=cos(thetamat).^elemfac;
elempatdB=20*log10(abs(elempat+eps));
%Compute Pattern
Pattern=AF.*elempat;
Patternmag=abs(Pattern);
Patternmagnorm=Patternmag/max(max(Patternmag)); %normalized pattern
PatterndB=20*log10(Patternmag+eps);
PatterndBnorm=20*log10(Patternmagnorm+eps);
%Set floor
dBfloor=-80;
PatterndB(find(PatterndB < dBfloor))=dBfloor;
PatterndBnorm(find(PatterndBnorm < dBfloor))=dBfloor;
%Generate Figures
figure 1)
clf
surf(u,v,PatterndB);hold on
plot(uo,vo,’yx’);
caxis([-50 20*log10(M*N)+5])
view(0,90),shading interp,colorbar
title([‘Pattern,’,num2str(M),’x’,num2str(N),’ Elements’]);
xlabel(‘u’)
ylabel(‘v’)
figure 11)
clf
surf(u,v,PatterndBnorm);hold on
plot(uo,vo,’yx’);
caxis([-80 0])
view(0,90),shading interp,colorbar
title([‘Pattern,’,num2str(M),’x’,num2str(N),’ Elements’]);
xlabel(‘u’)
ylabel(‘v’)
[row col]=find(PatterndB==max(max(PatterndB))); %find row and column of max
figure 2);clf;
plot(u(row,:),PatterndB(row,:),’b-’);
grid on;
title(‘Azimuth Cut of Array Pattern’);
xlabel(‘u’);
ylabel(‘dB’);
axis([-1 1 -20 80]);
figure 3);clf;
plot(v(:,col),PatterndB(:,col),’b-’);
grid on;
title(‘Elevation Cut of Array Pattern’);
xlabel(‘u’);
ylabel(‘dB’);
axis([-1 1 -20 80]);
figure 4)
clf
surf(u,v,elempatdB),hold
caxis([-20 0])
view(3),shading interp,colorbar
title([‘Element Pattern,’,num2str(M),’x’,num2str(N),’ Elements’])
plot(uo,vo,’yx’)
view(0,90);
xlabel(‘u’)
ylabel(‘v’)
figure 5)
clf
plot(x,y,’b.’);grid on;hold on;
switch removeElem
case 1
plot(x(DeleteList),y(DeleteList),’ro’)
case 2
[xplot yplot]=find(wgts==0);
plot(x(xplot,yplot),y(xplot, yplot),’ro’)
end
title(‘Locations of Array Elements’)
xlabel(‘X (in) ‘);
ylabel(‘ Y (in) ‘);
axis([-40 40 -40 40]);
return


