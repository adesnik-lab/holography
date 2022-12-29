function DE = DEfromSLMCoords(slmCoord)
load('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\ActiveCalib.mat')

SIcoord = function_SLMtoSI(slmCoord,CoC);
DE = computeDEfromList(SIcoord',{1},1);
