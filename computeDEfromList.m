function [AttenuationCoeffs, DElist] = computeDEfromList(SICoordinates,ROIs,weights)

%load current CoC
load('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\ActiveCalib.mat','CoC');

%catch targets that are below minimal diffraction efficiency 
DEfloor = 0.05;

[SLMCoordinates] = function_SItoSLM(SICoordinates',CoC)';
AttenuationCoeffs =SLMCoordinates(4,:);
lowDE = AttenuationCoeffs<DEfloor;
AttenuationCoeffs(lowDE)=DEfloor;
disp([num2str(sum(lowDE)) ' Target(s) below Diffraction Efficiency floor (' num2str(DEfloor) ').']);



%%

if size(weights,1)~=1
    weights=weights';
end

for i=1:numel(ROIs)
    ROIselection =ROIs{i};
    myattenuation = AttenuationCoeffs(ROIselection);
    energy = 1./myattenuation;
    energy = energy.*weights(ROIselection);
    energy = energy/sum(energy);
    DElist(i) = sum(energy.*myattenuation);
end