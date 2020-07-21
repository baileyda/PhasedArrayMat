%% Function to Compute 2D AF
% Arik D. Brown
function [AF, AF_mag, AF_dB, AF_dBnorm] =...
Compute_2D_AF(wgts,nelemsx,nelemsy,dx_in,dy_in,f_GHz,fo_GHz,...
randerror_amp,randerror_phs,u,v,uo,vo)
lambda=11.803/f_GHz;%wavelength(in)
lambdao=11.803/fo_GHz;%wavelength at tune freq(in)
k=2*pi/lambda;%rad/in
ko=2*pi/lambdao;%rad/in
AF=zeros(size(u));
phasex_mat=k*u-ko*uo;
phasey_mat=k*v-ko*vo;
randerror=randerror_amp.*exp(1j*randerror_phs*pi()/180);
for ii=1:nelemsx
dx=(ii-(nelemsx+1)/2)*dx_in;
for jj=1:nelemsy
dy=(jj-(nelemsy+1)/2)*dy_in;
AF = AF + wgts(ii,jj).*randerror(jj,ii)*...
exp(1j*dx*phasex_mat).*exp(1j*dy*phasey_mat);
end
end
[AF_mag AF_dB AF_dBnorm] = process_matrix(AF);
