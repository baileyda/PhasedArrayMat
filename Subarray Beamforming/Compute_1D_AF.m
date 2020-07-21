%% Function to Compute 1D AF
% Arik D. Brown
function [AF, AF_mag, AF_dB, AF_dBnorm] =...
Compute_1D_AF(wgts,nelems,d_in,f_GHz,fo_GHz,u,uo)
lambda=11.803/f_GHz;%wavelength(in)
lambdao=11.803/fo_GHz;%wavelength at tune freq(in)
k=2*pi/lambda;%rad/in
ko=2*pi/lambdao;%rad/in
AF=zeros(1,length(u));
for ii=1:nelems
AF = AF+wgts(ii)*exp(1j*(ii-(nelems+1)/2)*d_in*(k*u-ko*uo));
end
[AF_mag AF_dB AF_dBnorm] = process_vector(AF);