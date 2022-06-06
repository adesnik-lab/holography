
%% Extract data to fit OptotuneZ as a function of camera XYZ

disp('Fitting optotune to Camera... extracting optotune depths')
tFits=tic;

out.SIVals =SIVals;
out.XYSI =XYSI;
out.zsToUse =zsToUse;
out.SIUZ = SIUZ;

nGrids =size(SIVals,2);
nOpt = size(zsToUse,2);
fastWay = 0;

 clear SIpeakVal SIpeakDepth
fprintf('Extracting point: ')
parfor i=1:nGrids
    for k=1:nOpt
        if fastWay
            [a, b] = max(SIVals(:,i,k));
            SIpeakVal(i,k)=a;
            SIpeakDepth(i,k) =SIUZ(b);
        else
            try
                ff = fit(SIUZ', SIVals(:,i,k), 'gauss1');
                SIpeakVal(i,k) =ff.a1;
                SIpeakDepth(i,k) =ff.b1;
            catch
                SIpeakVal(i,k) = nan;
                SIpeakDepth(i,k) = nan;
            end
        end
    end
    fprintf([num2str(i) ' '])
    if mod(i,25)==0
        disp(' ')
    end
end

fprintf('\ndone\n')


b1 = SIpeakVal;
b2 = SIpeakDepth;

excl = SIpeakVal<(threshHold/1.5);
disp([num2str(sum(excl(:))) ' points excluded b/c below threshold'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;
% excl = SIpeakVal>20;
% SIpeakVal(excl)=nan;
% SIpeakDepth(excl)=nan;

excl = SIpeakDepth<-50;
disp([num2str(sum(excl(:))) ' points excluded b/c too deep'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;


%% CamToOpt
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
    2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

camXYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
camXYZ(3,:) =  SIpeakDepth(:);

camPower = SIpeakVal(:);

optZ = repmat(zsToUse,[nGrids 1]);
optZ = optZ(:);

testSet = randperm(numel(optZ),50);

otherSet = ones([numel(optZ) 1]);
otherSet(testSet)=0;
otherSet = logical(otherSet);

refAsk = (camXYZ(1:3,otherSet))';
refGet = optZ(otherSet);

camToOpto =  polyfitn(refAsk,refGet,modelterms);


Ask = camXYZ(1:3,testSet)';
True = optZ(testSet);

Get = polyvaln(camToOpto,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);

CoC.camToOpto= camToOpto;
%%fig
figure(201);clf
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],camPower,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Fluorescence intensity by space')
c = colorbar;
c.Label.String = 'Fluorescent Intensity';
axis square

figure(2);clf
subplot(2,2,1)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],optZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,2)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],polyvaln(camToOpto,camXYZ(1:3,:)'),'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Estimated Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,3)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],polyvaln(camToOpto,camXYZ(1:3,:)')-optZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Error (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,4)
c = sqrt((polyvaln(camToOpto,camXYZ(1:3,:)')-optZ).^2);
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],c,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Error RMS (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

%%optZtoCam
cam2XYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
cam2XYZ(3,:) =  optZ(:);
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

CoC.OptZToCam= OptZToCam;

out.CoC=CoC;
out.SIfitModelTerms = modelterms;
%%fig
figure(3);clf
subplot(2,2,1)
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],obsZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Measured Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,2)
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],polyvaln(OptZToCam,cam2XYZ(1:3,:)'),'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Estimated Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,3)
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],polyvaln(OptZToCam,cam2XYZ(1:3,:)')-obsZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Error (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,4)
c = sqrt((polyvaln(OptZToCam,cam2XYZ(1:3,:)')-obsZ).^2);
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],c,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Error RMS (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

disp(['All fits took ' num2str(toc(tFits)) 's']);
fitsT = toc(tFits);