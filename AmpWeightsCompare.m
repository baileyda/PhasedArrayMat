%% Plot Different Amplitude Weights
% Arik D. Brown
%% Enter Inputs
wgts.N=50;
wgts.nbar=5;
wgts.SLL_vec=[20 30 35];
wgts.uni_vec=ones(1,wgts.N);
wgts.tay_vec=[Taylor(wgts.N,wgts.SLL_vec(1),wgts.nbar)...
Taylor(wgts.N,wgts.SLL_vec(2),wgts.nbar)...
Taylor(wgts.N,wgts.SLL_vec(3),wgts.nbar)];
%% Plot Weights

figure 1),clf
plot(1:wgts.N,wgts.uni_vec,’-o’,’linewidth’,2,’color’,[0 0 1]),hold
plot(1:wgts.N,wgts.tay_vec(:,1),’--’,’linewidth’,2.5,’color’,[0 .7 0])
plot(1:wgts.N,wgts.tay_vec(:,2),’-.’,’linewidth’,2.5,’color’,[1 0 0])
plot(1:wgts.N,wgts.tay_vec(:,3),’:’,’linewidth’,2.5,’color’,[.7 0 1])
grid
% xlabel(‘Element Number’,’fontweight’,’bold’,’fontsize’,14)
% ylabel(‘Voltage’,’fontweight’,’bold’,’fontsize’,14)
set(gca,’fontweight’,’bold’,’fontsize’,14)
legend(‘Uniform Distribution’,’25 dB Taylor Distribution’,...
‘30 dB Taylor Distribution’,’35 dB Taylor Distribution’)
