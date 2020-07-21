%% Function to Compute 2D AF with Quantization
% Arik D. Brown
function [AF, AF_mag, AF_dB, AF_dBnorm] =...
Compute_2D_AF_quant(wgts,nelemsx,nelemsy,dx_in,dy_in,f_GHz,fo_GHz,...
randerror_amp,randerror_phs,nbits,u,v,uo,vo)
lambda=11.803/f_GHz;%wavelength(in)
lambdao=11.803/fo_GHz;%wavelength at tune freq(in)
k=2*pi/lambda;%rad/in
ko=2*pi/lambdao;%rad/in
AF=zeros(size(u));
xpos_vec=([1:nelemsy]-(nelemsx+1)/2)*dx_in;
ypos_vec=([1:nelemsy]-(nelemsy+1)/2)*dy_in;
[xpos_mat,ypos_mat]=meshgrid(xpos_vec,ypos_vec);
LSB=360/(2^nbits);
randerror=randerror_amp.*exp(1j*randerror_phs*pi()/180);
for ii=1:nelemsx
for jj=1:nelemsy
phase1=k*(xpos_mat(ii,jj)*u+ypos_mat(ii,jj)*v);
phase2=-ko*(180/pi)*(xpos_mat(ii,jj)*uo+ypos_mat(ii,jj)*vo);
phase2_quant1=phase2/LSB;
quant_delta=phase2_quant1-floor(phase2_quant1);
if quant_delta <= 0.5
phase2_quant2=floor(phase2_quant1)*LSB;
elseif quant_delta > 0.5
phase2_quant2=ceil(phase2_quant1)*LSB;
end
AF = AF + wgts(ii,jj).*randerror(jj,ii)*...
exp(1j*phase1).*exp(1j*phase2_quant2*pi/180);
end
end
[AF_mag AF_dB AF_dBnorm] = process_matrix(AF);
