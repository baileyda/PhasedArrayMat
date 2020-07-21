%
%% This Code Plots Beamwidth vs. Frequency and Scan Angle
% Arik D. Brown
%% Input Parameters
BW.k=0.886;%Beamwidth Factor (radians)
BW.f_vec=[1 5 10 15];%Frequency in GHZ
BW.lambda_vec=0.3./BW.f_vec;%meters
BW.L=1;%Aperture Length in meters
BW.thetao_vec=0:5:60;%Degrees
%% Calculate Beamwidths
[BW.lambda_mat BW.thetao_mat]=meshgrid(BW.lambda_vec,BW.thetao_vec);
BW.mat_rad=BW.k*BW.lambda_mat./(BW.L*cosd(BW.thetao_mat));
BW.mat_deg=BW.mat_rad*180/pi;
%% Plot
figure(1),clf
plot(BW.thetao_mat,BW.mat_deg,’linewidth’,2)
grid
set(gca,’fontsize’,16,’fontweight’,’b’)
xlabel(‘Scan Angle (Degrees)’,’fontsize’,16,’fontweight’,’b’)
ylabel(‘Beamwidth (degrees)’,’fontsize’,16,’fontweight’,’b’)
legend(‘1 GHz’,’5 GHz’,’10 GHz’,’15 GHz’)