%% Function to Compute Integrated Gain
% Arik D. Brown
function [PeakGain IdealGain PeakGaindB IdealGaindB GainPattern_mag...
GainPattern_dB GainPattern_dBnorm] =...
Compute_2D_INTGAIN(Pattern_mag,thetavec,phivec,thetamat,...
nelemx,nelemy,dx_in,dy_in,fGHz)
%% Compute Integrated Gain
numptstheta=length(thetavec);
numptsphi=length(phivec);
thetarad=thetavec*pi/180;
phirad=phivec*pi/180;
thetamat_rad=thetamat*pi/180;
dphi=(phirad(length(phirad))-phirad(1))/numptsphi;
dtheta=(thetarad(length(phirad))-thetarad(1))/numptstheta;
dsintheta=abs(sin(thetamat_rad));
GainPattern=(sum(sum((dsintheta*dphi*dtheta))))*2*Pattern_mag.^2/...
(sum(sum( (Pattern_mag.^2) .*(dsintheta*dphi*dtheta) )));
PeakGain=max(max(GainPattern));
PeakGaindB=10*log10(PeakGain);
[GainPattern_mag GainPattern_dB GainPattern_dBnorm] =...
process_matrix2(GainPattern);
%% Compute Ideal Gain
Area=nelemx*nelemy*dx_in*dy_in;
lambda=11.803/fGHz;
IdealGain=(4*pi*Area)/lambda.^2;
IdealGaindB=10*log10(IdealGain);
