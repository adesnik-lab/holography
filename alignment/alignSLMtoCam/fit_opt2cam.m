function result = fit_opt2cam(optotunePlanes, SIpeakDepth, nOpt, XYSI)

modelterms = optotuneModelTerms();

cam2XYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
cam2XYZ(3,:) =  optotunePlanes(:);
obsZ =  SIpeakDepth(:);

testSet = randperm(numel(obsZ),50);
otherSet = ones([numel(obsZ) 1]);
otherSet(testSet)=0;
otherSet = logical(otherSet);

refAsk = (cam2XYZ(1:3,otherSet))';
refGet = obsZ(otherSet);

OptZToCam =  polyfitn(refAsk,refGet,modelterms);


Ask = cam2XYZ(1:3,testSet)';
True = obsZ(testSet);

Get = polyvaln(OptZToCam,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp(['The mean error in Optotune depth prediction is : ' num2str(meanRMS) 'um']);
disp(['The Max error is: ' num2str(max(RMS)) 'um'])

result.OptZToCam = OptZToCam;
result.cam2XYZ = cam2XYZ;
result.obsZ = obsZ;
