function [SLMXYZP] = function_SItoSLM(SIXYZ,CoC)
%Function takes coordinates in SI space XY and optotuneZ (n x 3 matrix)
%and CoC as created in alignSLMtoCam
%Outputs in SLM coordinates XYZ and power normalization (n x 4 Matrix)
approachMethod=1;

SItoCam = CoC.SItoCam;
OptZToCam = CoC.OptZToCam;
camToSLM = CoC.camToSLM;
SLMtoPower = CoC.SLMtoPower;

estCamXY = function_Eval3DCoC(SItoCam,SIXYZ);
estCamZ = polyvaln(OptZToCam,estCamXY);
estCamXYZ = [estCamXY(:,1:2), estCamZ];
estSLM = function_Eval3DCoC(camToSLM,estCamXYZ);
estPower = polyvaln(SLMtoPower,estSLM);


if approachMethod ==2
%shortcut approach
estSLM = function_Eval3DCoC(CoC.SItoSLM,SIXYZ);
% estSI = function_Eval3DCoC(CoC.SLMtoSI,SLMXYZ);
end

SLMXYZP = [estSLM estPower];
