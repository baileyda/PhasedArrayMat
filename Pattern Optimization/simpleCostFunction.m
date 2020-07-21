function [Score]=simplePS_CostFunction(ComplexWgts,ElementWgtVoltage,LowerEnvelope,UpperEnvelope,AutomaticallyExemptMainBeam,EnforceSymmetry);
Nelements=prod(size(ComplexWgts));
if EnforceSymmetry
ComplexWgts(Nelements/2+1:Nelements)=flipud(ComplexWgts(1:Nelements/2));
end
NPatternPoints=prod(size(ElementWgtVoltage));
Pattern=abs(fftshift ( fft(ComplexWgts,NPatternPoints))).*ElementWgtVoltage;
PatternNorm=Pattern/max(Pattern);
LowerExcess=max(LowerEnvelope-PatternNorm,0);
UpperExcess=max(PatternNorm-UpperEnvelope,0);
if AutomaticallyExemptMainBeam
MaxValue=max(PatternNorm);
IndexOfMaxValue=find(PatternNorm==MaxValue);
Spikes=diff(sign(diff(PatternNorm)));
NullIndices=find(Spikes==2)+1;
PeakIndices=find(Spikes==-2)+1;
RightNullIndices=NullIndices(find(NullIndices>IndexOfMaxValue(1)));
LeftNullIndices=NullIndices(find(NullIndices<IndexOfMaxValue(1)));
if isempty(RightNullIndices)
IndexOfRightNull=NPatternPoints;
else
IndexOfRightNull=RightNullIndices(1);
end

if isempty(LeftNullIndices)
IndexOfLeftNull=1;
else
IndexOfLeftNull=LeftNullIndices(prod(size(LeftNullIndices)));
end
LowerExcess(IndexOfLeftNull:IndexOfRightNull)=0;
UpperExcess(IndexOfLeftNull:IndexOfRightNull)=0;
end
Score=10*log10(sum(LowerExcess)+sum(UpperExcess)+eps);

