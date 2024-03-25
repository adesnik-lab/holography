function result = fit_cam2opt(optotunePlanes, SIpeakDepth, SIpeakVal, XYSI)

ntest = 50;

%XY spatial calibration model for Power interpolations
modelterms = optotuneModelTerms();

nOpt = numel(optotunePlanes);
camXYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
camXYZ(3,:) =  SIpeakDepth(:);

camPower = SIpeakVal(:);
nGrids =size(SIVals,2);
optZ = repmat(optotunePlanes,[nGrids 1]);
optZ = optZ(:);

% generate test data
testSet = randperm(numel(optZ), ntest);

% select train data
otherSet = ones([numel(optZ) 1]);
otherSet(testSet) = 0;
otherSet = logical(otherSet);

% refAsk is cam position, XYZ
refAsk = (camXYZ(1:3,otherSet))';
% refGet is extracted opto position
refGet = optZ(otherSet);

% fit model
camToOpto =  polyfitn(refAsk, refGet, modelterms);
% Ask is requested, True is result
Ask = camXYZ(1:3,testSet)';
True = optZ(testSet);

% eval error
Get = polyvaln(camToOpto,Ask);
RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);


% results
result.camToOpto = camToOpto;
result.camXYZ = camXYZ;
result.camPower = camPower;
result.RMS = RMS;
result.meanRMS = meanRMS;


