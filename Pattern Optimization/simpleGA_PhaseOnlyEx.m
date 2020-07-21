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
% Taylor.m is included with the code in Chapters 1, 2, and 3
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
Nmarriages=25;
TournamentEligible=10;
TournamentCompetitors=5;
Ncrossovers=2;
MutationProbability=0.04;
MutationRangeMax=0.2;
MutationRangeDecay=1.5;
Norder=2*Nelements;
Nmutations=round(MutationProbability*Norder);
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
% Initialize Population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CurrentGeneration=rand(Norder,Npopulation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Score Population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:Npopulation
MagWgts=CurrentGeneration(1:Nelements,kk).*MagRange+MagMin;
PhsWgts=CurrentGeneration(Nelements+1:2*Nelements,kk).*PhsRange+PhsMin;
ComplexWgts=(MagWgts.*exp(i*PhsWgts));
CurrentGenerationScores(kk)=simpleCostFunction(...
ComplexWgts,ElementWgtVoltage,LowerEnvelope,UpperEnvelope,...
AutomaticallyExemptMainBeam,EnforceSymmetry);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Ready to Optimize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SortFitness,SortIndx]=sort(CurrentGenerationScores);
CurrentGeneration=CurrentGeneration(:,SortIndx);
CurrentGenerationScores=CurrentGenerationScores(SortIndx);
GlobalBestScore=CurrentGenerationScores(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Optimization Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:Niterations
t=(ii-1)/(Niterations-1);
MutationRange=MutationRangeMax*(exp(-MutationRangeDecay*t)-exp(-MutationRangeDecay));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Each Marriage:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:Nmarriages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose Mom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PermutedIndx=randperm(TournamentEligible);
CompetitorFitness=CurrentGeneration(PermutedIndx(1: TournamentCompetitors));
WinnerIndx=find(CompetitorFitness==min(CompetitorFitness));
Mom=CurrentGeneration(:,PermutedIndx(WinnerIndx(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose Dad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PermutedIndx=randperm(TournamentEligible);
CompetitorFitness=CurrentGeneration(PermutedIndx(1: TournamentCompetitors));
WinnerIndx=find(CompetitorFitness==min(CompetitorFitness));
Dad=CurrentGeneration(:,PermutedIndx(WinnerIndx(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Mom and Dad to Make Kids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Son=Dad;
Daughter=Mom;
PermutedIndx=randperm(Norder-1)+1;
CrossoverPoints=sort(PermutedIndx(1:Ncrossovers));
for CrossOverCount=1:Ncrossovers
Transfer=Son(CrossoverPoints(CrossOverCount):Norder);
Son(CrossoverPoints(CrossOverCount):Norder)=Daughter(CrossoverPoints(CrossOverCount):Norder);
Daughter(CrossoverPoints(CrossOverCount):Norder)=Transfer;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutate Son
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PermutedIndx=randperm(Norder);
MutationDecision=zeros(Norder,1);
MutationDecision(PermutedIndx(1:Nmutations))=1;
MutatedMax=min(Son+MutationRange/2, 1);
MutatedMin=max(Son-MutationRange/2, 0);
Mutations=rand(Norder,1).*(MutatedMax-MutatedMin)+MutatedMin;
Children(:,2*kk-1)=Son.*(1-MutationDecision)+Mutations.*MutationDecision;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutate Daughter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PermutedIndx=randperm(Norder);
MutationDecision=zeros(Norder,1);
MutationDecision(PermutedIndx(1:Nmutations))=1;
MutatedMax=min(Daughter+MutationRange/2, 1);
MutatedMin=max(Daughter-MutationRange/2, 0);
Mutations=rand(Norder,1).*(MutatedMax-MutatedMin)+MutatedMin;
Children(:,2*kk)=Daughter.*(1-MutationDecision)+Mutations.*MutationDecision;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Score Children
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=1:2*Nmarriages
MagWgts=Children(1:Nelements,kk).*MagRange+MagMin;
PhsWgts=Children(Nelements+1:2*Nelements,kk).*PhsRange+PhsMin;
ComplexWgts=(MagWgts.*exp(i*PhsWgts));
ChildrenScores(kk)=simpleCostFunction(...
ComplexWgts,ElementWgtVoltage,LowerEnvelope,UpperEnvelope,...
AutomaticallyExemptMainBeam,EnforceSymmetry);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Survival of the Fittest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CombinedGenerations=[CurrentGeneration Children];
CombinedScores=[CurrentGenerationScores ChildrenScores];
[SortFitness,SortIndx]=sort(CombinedScores);
SurvivorsIndx=SortIndx(1:Npopulation);
CurrentGeneration=CombinedGenerations(:,SurvivorsIndx);
CurrentGenerationScores=CombinedScores(SurvivorsIndx);
GlobalBestScore=CurrentGenerationScores(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Performance Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SavePopulationAverageScore(ii)=mean(CurrentGenerationScores);
SavePopulationBestScores(ii)=GlobalBestScore;
end
GlobalBestIndividual=CurrentGeneration(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Best Pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPatternPointsFine=NPatternPoints*8;
SineThetaFine=[-1:2/NPatternPointsFine:1-2/NPatternPointsFine];
CosineThetaFine=sqrt(1-SineThetaFine.^2);
ElementWgtPowerFine=CosineThetaFine.^CosineElementFactorExponent;
ElementWgtVoltageFine=sqrt(ElementWgtPowerFine).’;
ComplexWgts=(GlobalBestIndividual(1:Nelements).*MagRange+MagMin).*...
exp(i*(GlobalBestIndividual(Nelements+1:2*Nelements).*PhsRange+PhsMin));
if EnforceSymmetry
ComplexWgts(Nelements/2+1:Nelements)=flipud(ComplexWgts(1:Nelements/2));
end
ComplexWgts=ComplexWgts/max(abs(ComplexWgts));
BestPattern=fftshift (fft(ComplexWgts,NPatternPointsFine)).*ElementWgtVoltageFine;
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
plot([1:Niterations],SavePopulationBestScores,[1:Niterations],SavePopulationAverageScore,’LineWidth’,2)
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
title(‘Aperture weights from Genetic Algorithm’)
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

