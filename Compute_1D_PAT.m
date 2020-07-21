%% Function to Compute 1D PAT
% Arik D. Brown
function [PAT, PAT_mag, PAT_dB, PAT_dBnorm] =...
Compute_1D_PAT(EP,AF)
PAT=zeros(size(AF));
PAT=EP.*AF;
[PAT_mag PAT_dB PAT_dBnorm] =...
process_vector(PAT);
