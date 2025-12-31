function DE = DEfromSLMCoords(slmCoord)
load('C:\Users\Holography\Desktop\calibs\ActiveCalib.mat');

SIcoord = function_SLMtoSI(slmCoord,CoC);
DE = computeDEfromList(SIcoord',{1},1);
