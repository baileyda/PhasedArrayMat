clc
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nelements=100;
CosineElementFactorExponent=1.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excitation Constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor is included with the code in Chapters 1, 2, and 3
MagMin(1:Nelements,1)=Taylor(Nelements,30,5);
MagMax(1:Nelements,1)=Taylor(Nelements,30,5);
PhsMin(1:Nelements,1)=0;
PhsMax(1:Nelements,1)=pi/2;
EnforceSymmetry=0; % 1 for symmetric distrib, 0 for asymmetric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pattern Goals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UpperSidelobeGoal_u= [ -1 -.03 .03 .4 .4 .5 .5 1 ];
UpperSidelobeGoal_dB=[-55 -30 -30 -39.5 -60 -60 -42.1 -55];
LowerSidelobeGoal_u= [ -1 1];
LowerSidelobeGoal_dB=[-80 -80];
AutomaticallyExemptMainBeam=1; % 1 to exempt main beam from goals
NPatternPoints=512;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Niterations=1000;
Npopulation=50;
Phi1=2;
Phi2=2;
W=.4;
VRMSmax=.3;
Norder=2*Nelements;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Miscellaneous Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MagRange=MagMax-MagMin;
PhsRange=PhsMax-PhsMin;
SineTheta=[-1:2/NPatternPoints:1-2/NPatternPoints];
CosineTheta=sqrt(1-SineTheta.^2);
LowerEnvelope_dB=interp1(LowerSidelobeGoal_u+[0:size(LowerSidelobeGoal_u,2)-1]*eps, LowerSidelobeGoal_dB, SineTheta);
UpperEnvelope_dB=interp1(UpperSidelobeGoal_u+[0:size(UpperSidelobeGoal_u,2)-1]*eps, UpperSidelobeGoal_dB, SineTheta);
LowerEnvelope=10.^(LowerEnvelope_dB/20)’;
UpperEnvelope=10.^(UpperEnvelope_dB/20)’;
ElementWgtPower=CosineTheta.^CosineElementFactorExponent;
ElementWgtVoltage=sqrt(ElementWgtPower)’;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Particle Swarm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Swarm=rand(Norder,Npopulation);
SwarmVelocity=rand(Norder,Npopulation)*2-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Score Particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:Npopulation
MagWgts=Swarm(1:Nelements,kk).*MagRange+MagMin;
PhsWgts=Swarm(Nelements+1:2*Nelements,kk).*PhsRange+PhsMin;
ComplexWgts=(MagWgts.*exp(i*PhsWgts));
SwarmScores(kk)=simpleCostFunction(...
ComplexWgts,ElementWgtVoltage,LowerEnvelope,UpperEnvelope,...
AutomaticallyExemptMainBeam,EnforceSymmetry);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Ready to Optimize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SwarmScoreMemory=SwarmScores;
SwarmLocalBest=Swarm;
GlobalBestScore=min(SwarmScores);
IndexOfBest=find(SwarmScores==GlobalBestScore);
SwarmGlobalBest=Swarm(:,IndexOfBest(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Optimization Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:Niterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Each Particle:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:Npopulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SwarmSelfKnowledge=Phi1*rand(Norder,1).*...
(SwarmLocalBest(:,kk)-Swarm(:,kk));
SwarmSocialKnowledge=Phi2*rand(Norder,1).*...
(SwarmGlobalBest-Swarm(:,kk));
SwarmVelocity(:,kk)=W*SwarmVelocity(:,kk)+...
SwarmSelfKnowledge+SwarmSocialKnowledge;
VRMS=sqrt(sum(SwarmVelocity(:,kk).^2));
if VRMS>VRMSmax
SwarmVelocity(:,kk)=SwarmVelocity(:,kk)/VRMS*VRMSmax;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Swarm(:,kk)=Swarm(:,kk)+SwarmVelocity(:,kk);
Swarm(:,kk)=max(min(Swarm(:,kk),1),0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Score
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MagWgts=Swarm(1:Nelements,kk).*MagRange+MagMin;
PhsWgts=Swarm(Nelements+1:2*Nelements,kk).*PhsRange+PhsMin;
ComplexWgts=(MagWgts.*exp(i*PhsWgts));
SwarmScores(kk)=simpleCostFunction(...
ComplexWgts,ElementWgtVoltage,LowerEnvelope,UpperEnvelope,...
AutomaticallyExemptMainBeam,EnforceSymmetry);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Particle Memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SwarmScores(kk)<SwarmScoreMemory(kk)
SwarmScoreMemory(kk)=SwarmScores(kk);
SwarmLocalBest(:,kk)=Swarm(:,kk);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Group Knowledge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SwarmScores(kk)<GlobalBestScore
SwarmGlobalBest=Swarm(:,kk);
GlobalBestScore=SwarmScores(kk);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Performance Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SaveSwarmAverageScore(ii)=mean(SwarmScores);
SaveSwarmBestScores(ii)=GlobalBestScore;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Best Pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPatternPointsFine=NPatternPoints*8;
SineThetaFine=[-1:2/NPatternPointsFine:1-2/NPatternPointsFine];
CosineThetaFine=sqrt(1-SineThetaFine.^2);
ElementWgtPowerFine=CosineThetaFine.^CosineElementFactorExponent;
ElementWgtVoltageFine=sqrt(ElementWgtPowerFine).’;
ComplexWgts=(SwarmGlobalBest(1:Nelements).*MagRange+MagMin).*...
exp(i*(SwarmGlobalBest(Nelements+1:2*Nelements).*PhsRange+PhsMin));
if EnforceSymmetry
ComplexWgts(Nelements/2+1:Nelements)=flipud(ComplexWgts(1:Nelements/2));
end
ComplexWgts=ComplexWgts/max(abs(ComplexWgts));
BestPattern=fftshift (ff t(ComplexWgts,NPatternPointsFine)).*ElementWgtVoltageFine;
BestPattern_dB=20*log10(abs(BestPattern)+eps);
BestPatternNorm_dB=BestPattern_dB-max(BestPattern_dB);
ElapsedTimeMinutes=toc/60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
figscale=70; figoffsetx=20; figoffsety=20;
set(gcf,’Position’,[figoffsetx figoffsety round(11.25*figscale+figoffsetx) round(6.75*figscale+figoffsety)])
fontsize=12;
subplot(2,2,1)
set(gca,’FontSize’,fontsize)
plot([1:Niterations],SaveSwarmBestScores,[1:Niterations],SaveSwarmAverageScore,’LineWidth’,2)
xlabel(‘Iteration #’)
ylabel(‘Cost Function (arbitrary units)’)
legend(‘Best’,’Average’)
axis tight
subplot(2,2,2)
[ax,h1,h2]=plotyy([1:Nelements],abs(ComplexWgts),[1:Nelements],angle(ComplexWgts)-min(angle(ComplexWgts)));
set(h1,’Marker’,’.’,’LineStyle’,’none’);
set(h2,’Marker’,’.’,’LineStyle’,’none’);
axes(ax(1));
set(gca,’FontSize’,fontsize)
ylabel(‘Amplitude (Volts) (blue)’)
ylim([0 1])
axes(ax(2));
set(gca,’FontSize’,fontsize)
ylabel(‘Phase (radians) (green)’)
title(‘Aperture weights from Particle Swarm’)
subplot(2,1,2)
set(gca,’FontSize’,fontsize)
plot(SineThetaFine,BestPatternNorm_dB,’LineWidth’,2);
hold on
plot(LowerSidelobeGoal_u, LowerSidelobeGoal_dB, ‘r-.’,’LineWidth’,2)
plot(UpperSidelobeGoal_u, UpperSidelobeGoal_dB, ‘r:’,’LineWidth’,2)
axis tight
ylim([-80 0])
grid on
title([‘Far field Pattern (Target sidelobe constraints in red). Elapsed time = ‘ num2str(ElapsedTimeMinutes) ‘ minutes.’])
ylabel(‘Magnitude (dB)’)
xlabel(‘Sin(\theta)’)



