%% Function to Compute 1D EP
% Arik D. Brown
function [EP, EP_mag, EP_dB, EP_dBnorm] =...
Compute_1D_EP(theta_deg,EF)
EP=zeros(size(theta_deg));
EP=(cosd(theta_deg).^(EF/2));                      % Volts
[EP_mag, EP_dB, EP_dBnorm] = process_vector(EP);   
