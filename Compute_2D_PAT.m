%% Function to Compute 2D PAT
% Arik D. Brown
function [PAT, PAT_mag, PAT_dB, PAT_dBnorm] =...
Compute_2D_PAT(EP,AF)
PAT=zeros(size(AF));
PAT=EP.*AF;
[PAT_mag PAT_dB PAT_dBnorm] =...
process_matrix(PAT);
