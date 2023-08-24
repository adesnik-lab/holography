function alignSLMtoCam_ver_1_5()

%% Pathing and setup

rt = 'C:\Users\Holography\Desktop\holography\';
cd(rt)
addpath(genpath('C:\Users\Holography\Desktop\holography'))
addpath(genpath('C:\Users\Holography\Desktop\meadowlark'));
f5_idiot_check()

clear
close all
clc

tBegin = tic;

disp('Setting up stuff...');

Setup = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
Setup.useGPU =1;
Setup.SLM.is_onek = 1;

if Setup.useGPU
    disp('Getting gpu...'); %this can sometimes take a while at initialization
    g= gpuDevice;
end

% delete(gcp('nocreate'));
% parpool('IdleTimeout', 360);

slm = MeadowlarkSLM();
slm.start()

sutter = sutterController();
mvShortWait = 0.2;
mvLongWait = 2;

bas = bascam();
bas.start()
castImg = bas.castFun;
castAs = bas.castAs;

% reset rng to be more random
rng shuffle

disp('Ready! Good luck...')

%% look for objective in 1p

bas.preview()

%% connect to DAQ computer

% run this first then code on daq
% daq code is probably called alignCodeDAQ.m

laser = LaserSocket(42134);

%% connect to ScanImage computer

% run this code first, then 'autoCalibSI' on SI computer
fprintf('Waiting for msocket communication to ScanImage Computer... ')
% then wait for a handshake
srvsock2 = mslisten(42047);
SISocket = msaccept(srvsock2,15);
msclose(srvsock2);
sendVar = 'A';
mssend(SISocket, sendVar);

invar = [];
while ~strcmp(invar,'B')
    invar = msrecv(SISocket,.1);
end
fprintf('done.\r')

%% set power levels

% highly dependent on previous/current power calibration and alignment
% generally expect 10 mW for full power, 50 mW 10x divided mode, and 70 mW
% at 100x divided mode
%
% put a hologram near the zero-order and set power so that it is roughly
% 80% of camera max
%
% this power will be used throughout the calibration and is appropriately
% scaled for multi-target holograms and hole-burning

pwr = 20;
slmCoords = [.4 .4 0.01 1];

disp(['Individual hologram power set to ' num2str(pwr) 'mW.'])

DEestimate = DEfromSLMCoords(slmCoords);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);

slm.feed(Holo)
laser.set_power(pwr)
bas.preview()
laser.set_power(0)

%% check pixel val
bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

laser.set_power(pwr)
data = bas.grab(10);
laser.set_power(0)

data = mean(data,3);
frame = max(data-bgd, 0);
[x,y] = function_findcenter(frame);

frame_crop = data(x-20:x+20, y-20:y+20);

figure(6)
clf
imagesc(frame_crop)
axis square
colorbar
% clim([0 255])
xL=xlim;
yL=ylim;
mx = max(frame_crop,[],"all");
str = ['Peak Intensity: ' num2str(mx) ' A.U.'];
text(0.03*xL(2),0.03*yL(2),str,'HorizontalAlignment','left','VerticalAlignment','top', 'Color','w')

%% center FOV

disp('Find Focal Plane, Center and Zero the Sutter')
disp('Leave the focus at Zoom 1. at a power that is less likely to bleach (14% 25mW)') %25% 8/16/19
disp('Don''t forget to use Ultrasound Gel on the objective so it doesn''t evaporate')
laser.set_power(0)

bas.preview()

input('Turn off Focus and press any key to continue');

sutter.setRef()

mssend(SISocket,[0 0]);

disp('Make Sure the DAQ computer is running testMultiTargetsDAQ. and the SI computer running autoCalibSI');
disp('also make those names better someday')
disp('Make sure both lasers are on and the shutters open')
disp('Scanimage should be idle, nearly in plane with focus. and with the gain set high enough to see most of the FOV without saturating')

sutter.moveZ(100)

disp('testing the sutter double check that it moved to reference +100');
input('Ready to go (Press any key to continue)');

sutter.moveToRef()

tManual = toc(tBegin);

%% Create initial holos for coarse search

npts = 250;

slmXrange = [0.05 0.98];
slmYrange = [0.05 0.98];
slmZrange = [-0.015 0.015];

slmCoordsInitial = generateRandomHolos(slmXrange, slmYrange, slmZrange, npts);

newfig('Initial Coarse Holograms')
scatter3(slmCoordsInitial(1,:), ...
    slmCoordsInitial(2,:), ...
    slmCoordsInitial(3,:))
title('Initial Holograms')
drawnow

fprintf('Compiling initial holograms... ')
hololist = compileSingleTargetHolos(Setup, slmCoordsInitial);
fprintf('done.\r')

%% Measure background noise in camera

nBackgroundFrames = 10;

bgd_frames = bas.grab(nBackgroundFrames);
bgd = mean(bgd_frames, 3);

meanBgd = mean(bgd_frames, 'all');
stdBgd = std(single(bgd_frames), [], 'all');

newfig('Camera background')
imagesc(bgd)


%% Scan Image Planes Calibration
disp('Begining SI Depth calibration, we do this first in case spots burn holes with holograms')
clear SIdepthData

zsToUse = linspace(0,70,15); % ScanImage/Optotune planes
SIUZ = -15:5:130; % sutter planes
framesToAcquire = 10; % camera average
nframesCapture = framesToAcquire;

SIpts = numel(SIUZ);

%generate xy grid
%this is used bc SI depths are not necessarily parallel to camera depths
%but SI just generates a sheet of illumination
%therefore, we want to check a grid of spots that we will get depth info
%for
sz = size(bgd);
gridpts = 25;
xs = round(linspace(1,sz(1),gridpts+2));
ys = round(linspace(1,sz(2),gridpts+2));
xs([1 end])=[];
ys([1 end])=[];

clear dimx dimy XYSI
range = 15;
c=0;
for i=1:gridpts
    for k=1:gridpts
        c=c+1;
        dimx(:,c) = xs(i)-range:xs(i)+range;
        dimy(:,c) = ys(k)-range:ys(k)+range;
        XYSI(:,c) = [xs(i) ys(k)];
    end
end

disp(['We will collect ' num2str(numel(zsToUse)) ' planes.'])

% for each optotune plane, locate in basler with sutter Z

tSI = tic;

SIVals = zeros([SIpts c numel(zsToUse)]);

for k =1:numel(zsToUse)
    t = tic;

    z = zsToUse(k);

    disp(['Testing plane ' num2str(z) ': ']);

    mssend(SISocket,[z 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(SISocket,0.01);
    end

    dataUZ = zeros([sz SIpts]);
    for i = 1:numel(SIUZ)
        fprintf([num2str(round(SIUZ(i))) ' ']);

        sutter.moveZ(SIUZ(i))

        if i==1
            pause(mvLongWait)
        else
            pause(mvShortWait);
        end

        data = bas.grab(nframesCapture);

        % post process capture
        data = mean(data, 3);
        frame = max(data-bgd, 0);
        frame = imgaussfilt(frame, 2);

        dataUZ(:,:,i) =  frame;

    end
    sutter.moveToRef()
    pause(mvLongWait)

    mssend(SISocket,[z 0]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(SISocket,0.01);
    end

    figure(1212);
    for hbtmpi=1:numel(SIUZ)
        subplot(6,6,hbtmpi); colorbar;
        imagesc(dataUZ(:,:,hbtmpi));
    end
    pause(.5)

    for i =1:c
        SIVals(:,i,k) = squeeze(mean(mean(dataUZ(dimx(:,i),dimy(:,i),:))));
    end

    disp([' Took ' num2str(toc(t)) 's']);
end

mssend(SISocket,'end');

disp(['Scanimage calibration done whole thing took ' num2str(toc(tSI)) 's']);
siT=toc(tSI);

%%
tFits=tic;
disp('No Longer Putting off the actual analysis until later, Just saving for now')
out.SIVals =SIVals;
out.XYSI =XYSI;
out.zsToUse =zsToUse;
out.SIUZ = SIUZ;

save('TempSIAlign.mat','out')

%% First Fits
%Extract data to fit OptotuneZ as a function of camera XYZ

disp('Fitting optotune to Camera... extracting optotune depths')


out.SIVals =SIVals;
out.XYSI =XYSI;
out.zsToUse =zsToUse;
out.SIUZ = SIUZ;
nGrids =size(SIVals,2);
nOpt = size(zsToUse,2);
fastWay = 01;

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
%%
SIpeakVal = b1;
SIpeakDepth = b2;

SIThreshHoldmodifier = 1.5;
SIThreshHold =SIThreshHoldmodifier*stdBgd/sqrt(nBackgroundFrames + framesToAcquire);

excl = SIpeakVal < SIThreshHold;
disp([num2str(numel(SIpeakDepth)) ' points total before exclusions'])
disp([num2str(sum(excl(:))) ' points excluded b/c below threshold'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;

excl = SIpeakVal>255;
excl = SIpeakDepth<-10 | SIpeakDepth > 130; %upper bound added 7/15/2020 -Ian

disp([num2str(sum(excl(:))) ' points excluded b/c too deep'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;
disp([num2str(sum(~isnan(SIpeakDepth), 'all')) ' points remaining'])


%%CamToOpt
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
f201=figure(201);clf
f201.Units = 'normalized';
f201.Position = [0 0 0.25 0.45];
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],camPower,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Fluorescence intensity by space')
c = colorbar;
c.Label.String = 'Fluorescent Intensity';
axis square

f2 = figure(2);clf
f2.Units = 'normalized';
f2.Position = [0 0.4 0.4 0.4];
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
f3 = figure(3);clf
f3.Units = 'normalized';
f3.Position = [0.4 0.4 0.4 0.4];
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


%% Coarse Data
tstart=tic;

disp('Begining Coarse Holo spot finding')

numFramesCoarseHolo = 5; % number frames to collect
coarsePts = 9; % odd number please
coarseUZ = linspace(-25,150,coarsePts);

laser.set_power(0)
laser.flush()

vals = nan(coarsePts,npts);
xyLoc = nan(2,npts);

sz = size(bgd);
sizeFactor = 1;
newSize = sz / sizeFactor;

dataUZ2 = zeros([newSize  numel(coarseUZ) npts], castAs);
maxProjections = castImg(zeros([newSize  npts]));

for i = 1:numel(coarseUZ)
    fprintf(['First Pass Holo, Depth: ' num2str(coarseUZ(i)) '. Holo : '])
    t = tic;

    sutter.moveZ(coarseUZ(i))

    if i==1
        pause(mvLongWait)
    else
        pause(mvShortWait)
    end

    for k=1:npts
        fprintf([num2str(k) ' ']);

        if mod(k,25)==0
            fprintf('\n')
        end

        % send holo
        slm.feed(hololist(:,:,k))

        % show holo and grab images
        laser.set_power(pwr)
        frame = bas.grab(numFramesCoarseHolo);
        laser.set_power(0)

        % post process
        frame = mean(frame, 3);
        frame = max(frame-bgd, 0);
        frame = imgaussfilt(frame, 2);
        frame = imresize(frame, newSize);
        dataUZ2(:,:,i,k) =  frame;

    end
    fprintf(['\nPlane Took ' num2str(toc(t)) ' seconds\n'])

end

% cast to unit8
dataUZ2 = uint8(dataUZ2);

sutter.moveToRef()

%%
disp('Calculating Depths and Vals')

range = ceil(5/sizeFactor);

for k=1:npts
    dataUZ = dataUZ2(:,:,:,k);
    mxProj = max(dataUZ,[],3);
    [ x,y ] = function_findcenter(mxProj );
    xyLoc(:,k) = [x,y]*sizeFactor;

    maxProjections(:,:,k) = mxProj;

    dimx = max((x-range),1):min((x+range), size(mxProj,1));
    dimy = max((y-range),1):min((y+range), size(mxProj,2));

    thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
    vals(:,k) = thisStack;
    depthIndex = find(thisStack == max(thisStack),1);

    fprintf(['Spot ' num2str(k) ' centered at depth ' num2str(round(coarseUZ(depthIndex)))...
        'um. Value: ' num2str(round(vals(depthIndex,k))) '\n']);
end

fprintf(['All Done. Total Took ' num2str(toc(tstart)) 's\n']);

%% Second pass, multi-target version
% this section takes the holograms X,Y,Z coarse locations from the first
% step and then performs a fine search for them based off of thier
% approximate z plane, only done on the initialized points

% assign search params
finePts = 13; % odd number please
fineRange = 40;
coaseIncludeThreshScalar = 3;
nframesCapture = 10;

% holo generation params
iterationsBeforeStop = 1000;
distanceThreshold = 30; % spacing between points, changed from 50 on 7/15/20 bc new cam
size_of_holo = 10; % number of points in holo
cores = 2; % for parallel compute of holos

disp('Begin multi-target z search.')
multi_time = tic;

laser.flush()

% Generate multi-target holos based off coarse search data
targ_time = tic;

[coarseVal, coarseZidx] =max(vals,[],1);
zDepthVal = coarseUZ(coarseZidx);
zdepths = unique(zDepthVal);
n_planes = numel(zdepths);

% inclusion threshold added based on frames acquired; more stringent then SI.
coarseInclusionThreshold = coaseIncludeThreshScalar * stdBgd/sqrt(numFramesCoarseHolo + nBackgroundFrames);
zDepthVal(coarseVal < coarseInclusionThreshold) = NaN;

xyzLoc = [xyLoc; zDepthVal]; %fix this later (?)

clear slmMultiCoords basCoords targ_list targListIndiv slmMultiCoordsIndiv tempTargList
disp('Making holos.')
for i=1:n_planes % this will be the number of holograms
    % index the z depth
    z = zdepths(i);
    targ_idx = find(xyzLoc(3,:)==z);
    slmMultiCoords{i} = slmCoordsInitial(:,targ_idx);
    basCoords{i} = xyzLoc(:,targ_idx);
    targ_list{i} = targ_idx;
end

for i=1:n_planes
    % get real and slm coords from coarse
    dist = pdist2(basCoords{i}',basCoords{i}');
    temp = rand(size(dist,1));
    dist(find(diag(diag(temp)))) = nan;  %#ok<FNDSB>
    tempTargList = 1:numel(targ_list{i});
    iterCounter = 0;
    multiHoloCounter = 0;
    keepGoing = 1;
    doThisOnce = 0;
    slmMultiCoordsIndiv{i} =[];
    targListIndiv{i}=[];

    while keepGoing
        iterCounter=iterCounter+1;
        if numel(tempTargList) <= size_of_holo
            testIdx = tempTargList;
            IdxofTempTargetList = 1:numel(tempTargList);
            keepGoing = 0;
        else
            IdxofTempTargetList = randperm(numel(tempTargList), size_of_holo);
            testIdx = tempTargList(IdxofTempTargetList);
        end

        %test if good
        subDist = dist(testIdx, testIdx);
        if any(subDist(:) < distanceThreshold)
            good =0;
        else
            good = 1;
        end

        if good
            multiHoloCounter = multiHoloCounter + 1;
            slmMultiCoordsIndiv{i}{multiHoloCounter} = slmMultiCoords{i}(:, testIdx);
            targListIndiv{i}{multiHoloCounter} = targ_list{i}(testIdx);
            iterCounter = 0;
            tempTargList(IdxofTempTargetList)=[];
        else
            if iterCounter > iterationsBeforeStop && doThisOnce
                keepGoing=0;
            elseif iterCounter > iterationsBeforeStop
                size_of_holo=max(round(size_of_holo/2), 3);
                iterCounter = 0;
                doThisOnce = 1;
            end
        end
    end
end

%% make the holos
disp('Setting up compute for multi-targets...');
Setup.CGHMethod = 2;
Setup.GSoffset = 0;
Setup.verbose = 0;
Setup.useGPU = 1;

if cores > 1
    p = gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=cores
        delete(p);
        parpool(cores);
    end
end

clear slmShootCoords
holo_time = tic;
disp('Compiling holograms...')
planes = numel(slmMultiCoordsIndiv);

for i=1:planes
    pt = tic;
    holos_this_plane = numel(slmMultiCoordsIndiv{i});

    for k=1:holos_this_plane
        [ mtholo, ~, ~ ] = function_Make_3D_SHOT_Holos(Setup,slmMultiCoordsIndiv{i}{k}');
        mtholo_temp(k,:,:) = mtholo;
    end

    multiHolos{i} = mtholo_temp;
    disp(['Plane ' num2str(i) ' of ' num2str(planes) ' done!  Took ' num2str(toc(pt)) 's'])

end
disp(['Done. Took ' num2str(toc(holo_time)) 's'])


disp(['Overall took ' num2str(toc(targ_time)) 's to compile multi target holos'])

out.hololist = hololist;
out.slmCoords = slmCoordsInitial;

%%
clear peakValue peakDepth peakFWHM

box_range = 20; % 7/15/20 changed from 50 to 20 distance threshold is set to 50, this must be less to avoid trying to fit 2 holos
disp('shootin!')

% for every  plane
for i = 1:planes
    plane_time = tic;
    holos_this_plane = numel(slmMultiCoordsIndiv{i});

    disp(['Plane ' num2str(i) ' of ' num2str(planes)])

    % for every holo on that plane
    for j = 1:holos_this_plane
        holo_time = tic;
        disp(['Multi-target holo ' num2str(j) ' of ' num2str(holos_this_plane)])

        if size(slmMultiCoordsIndiv{i}{j},2) == 0 || size(slmMultiCoordsIndiv{i}{j},2) == 2 % or <3 ??
            continue
        end

        multi_pwr = size(slmMultiCoordsIndiv{i}{j},2) * pwr;

        target_ref = targListIndiv{i}{j}(1);
        expected_z = xyzLoc(3, target_ref);
        %expected_xy = xyzLoc(1:2, target_ref);
        %expected_xyz = xyzLoc(:,target_ref);

        fineUZ = linspace(expected_z-fineRange, expected_z+fineRange, finePts);
        dataUZ = nan([size(bgd(:,:,1))  finePts]);

        slm.feed(multiHolos{i}(j,:,:))

        fprintf('Depth: ')
        % for every sutter z plane
        for k = 1:finePts
            fprintf([num2str(round(fineUZ(k))) ' ']);

            sutter.moveZ(fineUZ(k))

            if i==1
                pause(mvLongWait)
            else
                pause(mvShortWait)
            end

            laser.set_power(multi_pwr)
            data = bas.grab(nframesCapture);
            laser.set_power(0)

            % subtract the background and filter
            frame = mean(data, 3);
            frame = max(frame-bgd, 0);
            frame = imgaussfilt(frame, 2);

            % store into dataUZ(x,y,z-plane)
            dataUZ(:,:,k) =  frame;
            dataUZ3{i}{j} = uint8(dataUZ);
            fineUZ3{i}{j} = fineUZ;
        end

        sutter.moveToRef()

        % OK, now parse the basler data in expected holo spots
        for targ = 1:size(slmMultiCoordsIndiv{i}{j},2)

            target_ref = targListIndiv{i}{j}(targ);
            expected_xyz = xyzLoc(:, target_ref);
            [x, y] = size(bgd);

            targX = expected_xyz(1)-box_range:expected_xyz(1)+box_range;
            targY = expected_xyz(2)-box_range:expected_xyz(2)+box_range;

            if max(targX)>x
                targX = expected_xyz(1)-box_range:x;
            end
            if max(targY)>y
                targY = expected_xyz(2)-box_range:y;
            end
            if min(targX)<1
                targX = 1:expected_xyz(1)+box_range;
            end
            if min(targY)<1
                targY = 1:expected_xyz(2)+box_range;
            end

            % catch for small boxes!!!!
            try
                % method 1 - rely on XY from first step
                targ_stack = squeeze(max(max(dataUZ(targX,targY,:))));
                mxProj = max(dataUZ(targX,targY,:), [], 3);
                [holo_x, holo_y] = function_findcenter(mxProj);
                xyFine{i}{j}(:,targ) = [holo_x, holo_y];
            catch
                targ_stack = nan(finePts,1);
                xyFine{i}{j}(targ) = nan;
            end

            % fit data
            try
                ff = fit(fineUZ', targ_stack, 'gauss1');
                peakValue{i}{j}(targ) = ff.a1;
                peakDepth{i}{j}(targ) = ff.b1;
                peakFWHM{i}{j}(targ) = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
            catch
                disp(['Error on fit! Holo: ', num2str(j), ' Target: ', num2str(targ)])
                peakValue{i}{j}(targ) = NaN;
                peakDepth{i}{j}(targ) = NaN;
                peakFWHM{i}{j}(targ) = NaN;
            end
        end

        fprintf('\n')
        disp(['Holo ' num2str(j) ' took ' num2str(toc(holo_time)) 's'])
    end

    disp(['Plane ' num2str(i) ' took ' num2str(toc(plane_time)) 's'])
end

fprintf(['All Done. Total Took ' num2str(toc(multi_time)) 's\n']);

sutter.moveToRef()
pause(mvLongWait)

%%
%% reshape matrix of values into vectors
tIntermediateFine = tic;

peakValueList = peakValue(:);
peakDepthList = peakDepth(:);
peakFWHMList = peakFWHM(:);

% current implementation does not threshold, should be included but more
% complicated here
c=0;
clear slmXYZ basXYZ1 basVal1 FWHMval
for i=1:planes

    %changed to account for when it skipped asked holos bc they were too small -Ian 7/16/2020
    if numel(peakDepth)<i
        holos_this_plane = 0;
    else
        holos_this_plane = numel(peakDepth{i});
    end

    for j=1:holos_this_plane
        for targ = 1:size(slmMultiCoordsIndiv{i}{j},2)
            c = c+1;
            slmXYZ(:,c) = slmMultiCoordsIndiv{i}{j}(:,targ);
            target_ref = targListIndiv{i}{j}(targ);

            %             basXYZ1_fine(1:2,c) = xyFine{i}{j}(targ);
            basXYZ1(1:2,c) = xyzLoc(1:2,target_ref);
            basXYZ1(3,c) = peakDepth{i}{j}(targ);
            basVal1(c) = peakValue{i}{j}(targ);
            FWHMval(c) = peakFWHM{i}{j}(targ);
        end
    end
end

basXYZ = basXYZ1;
basVal = basVal1;

%%

slmXYZBackup2 = slmXYZ;
basXYZBackup2 = basXYZ;
basValBackup2 = basVal;
FWHMBackup2   = FWHMval; %added 7/20/2020 -Ian

%% --- BEGIN INTERMEDIATE FITS ---- %%

%% exclude trials

slmXYZ = slmXYZBackup2;
basXYZ = basXYZBackup2;
basVal = basValBackup2;
FWHMVal = FWHMBackup2;

newfig('Bas Data')
subplot(1,2,1)
scatter3(basXYZ(1,:), basXYZ(2,:), basXYZ(3,:))
subplot(1,2,2)
histogram(basVal)


excludeTrials = all(basXYZ(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

excludeTrials = excludeTrials | basVal > bas.camMax;

basDimensions = size(bgd);
excludeTrials = excludeTrials | basXYZ(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ(3,:)<-5;
excludeTrials = excludeTrials | basXYZ(3,:)>200;


excludeTrials = excludeTrials | any(isnan(basXYZ(:,:)));
excludeTrials = excludeTrials | basVal < 1;
excludeTrials = excludeTrials | basVal > 250;
% excludeTrials = excludeTrials | basVal>(mean(basVal)+3*std(basVal));

slmXYZBackup = slmXYZ(:,~excludeTrials);
basXYZBackup = basXYZ(:,~excludeTrials);
basValBackup = basVal(:,~excludeTrials);
FWHMValBackup = FWHMVal(~excludeTrials); % added 7/20/2020 -Ian
disp(['num total holos: ' num2str(size(excludeTrials,2))])
disp(['num exlc holos: ' num2str(sum(excludeTrials))])
%%
f41=figure(41);
clf(41)
f41.Units = 'Normalized';
f41.Position = [0.05 0.4 0.5 0.5];
subplot(1,2,1)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 50, 'k', 'Filled', 'MarkerFaceAlpha',0.5)
hold on
scatter3(xyzLoc(1,:), xyzLoc(2,:), xyzLoc(3,:), 70,  'r', 'Filled', 'MarkerFaceAlpha', 0.7)
legend('Fine','Coarse')
title('Detected basXYZs')
subplot(1,2,2)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 75, basValBackup, 'Filled')
colorbar
colormap default
title('basXYZ and basVals (fine)')

f42=figure(42);
clf(42)
f42.Units = 'Normalized';
f42.Position = [0.05 0 0.40 0.5];
hold on
basx = 1:size(basVal,2);
plot(basx, basVal,'ro')
plot(basx(~excludeTrials), basVal(~excludeTrials), 'o')
legend('Excluded','Included')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')

%% fit SLM to Camera
% use model terms

% params
holdback = 30;
errScalar = 2.5; % should be between 2 and 3
% pxPerMu = 0.57723; % i think this is actually muPerPx, functionally
pxPerMu = 1;

basXYZ = basXYZBackup;
slmXYZ = slmXYZBackup;
basVal = basValBackup;

disp('Fitting SLM to Camera')
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1 ; 1 1 1 ;...
    2 0 0; 0 2 0; 0 0 2;  ...
    2 0 1; 2 1 0; 0 2 1; 1 2 0; 0 1 2;  1 0 2; ... ];  %XY spatial calibration model for Power interpolations
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2;];

% shuffle
reOrder = randperm(size(slmXYZ,2));
slmXYZ = slmXYZ(:,reOrder);
basXYZ = basXYZ(:,reOrder);

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = (basXYZ(1:3,1:end-holdback))';

figure(1286)
clf
ax = gca();
[SLMtoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,ax);
title('SLM to Cam v1')

Ask = refAsk;
True = refGet;
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(103);clf
subplot(1,2,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
% legend('Measured targets', 'Estimated Targets');
title({'Reference Data'; 'SLM to Camera'})

refRMS = sqrt(sum((Get-True).^2,2));
subplot(1,2,2)
scatter3(True(:,1),True(:,2),True(:,3),[],refRMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title({'Reference Data'; 'RMS Error in position'})
caxis([0 30])


Ask = (slmXYZ(1:3,end-holdback:end))';
True = (basXYZ(1:3,end-holdback:end))';
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(101);clf
subplot(1,3,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')


ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
legend('Measured targets', 'Estimated Targets');
title('SLM to Camera')

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp('Error based on Holdback Data...')
disp(['The RMS error: ' num2str(meanRMS) ' pixels for SLM to Camera']);


disp(['Thats approx ' num2str(meanRMS/pxPerMu) ' um']);

xErr = sqrt(sum((Get(:,1)-True(:,1)).^2,2));
yErr = sqrt(sum((Get(:,2)-True(:,2)).^2,2));
zErr = sqrt(sum((Get(:,3)-True(:,3)).^2,2));

disp('Mean:')
disp(['X: ' num2str(mean(xErr)/pxPerMu) 'um. Y: ' num2str(mean(yErr)/pxPerMu) 'um. Z: ' num2str(mean(zErr)) 'um.']);
disp('Max:')
disp(['X: ' num2str(max(xErr)/pxPerMu) 'um. Y: ' num2str(max(yErr)/pxPerMu) 'um. Z: ' num2str(max(zErr)) 'um.']);

subplot(1,3,2)
scatter3(True(:,1),True(:,2),True(:,3),[],RMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title('RMS Error in position')

refAsk = (basXYZ(1:3,1:end-holdback))';
refGet = (slmXYZ(1:3,1:end-holdback))';

camToSLM = function_3DCoC(refAsk,refGet,modelterms);

figure(1401);
clf
ax = gca();
[camToSLM, ~] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,ax);
title('cam to slm v1')

Ask = (basXYZ(1:3,end-holdback:end))';
True = (slmXYZ(1:3,end-holdback:end))';
Get = function_Eval3DCoC(camToSLM,Ask);

figure(101);
subplot(1,3,3)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Depth units')
legend('Measured targets', 'Estimated Targets');
title('Camera to SLM')

CoC.camToSLM=camToSLM;
CoC.SLMtoCam = SLMtoCam;

out.CoC=CoC;
out.CoCmodelterms = modelterms;

rtXYZ = function_Eval3DCoC(SLMtoCam,function_Eval3DCoC(camToSLM,basXYZ(1:3,end-holdback:end)'));

err = sqrt(sum((rtXYZ - basXYZ(1:3,end-holdback:end)').^2,2));
meanRTerr = nanmean(err);
disp(['The Mean Round Trip RMS error: ' num2str(meanRTerr) ' pixels (' num2str(meanRTerr/pxPerMu) ' um) camera to SLM to camera']);

%% fit power as a function of SLM
disp('Fitting Power as a function of SLM')

slmXYZ = slmXYZBackup;
basVal = basValBackup;

modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...

intVal = basVal;
intVal = sqrt(intVal); % convert fluorescence intensity (2P) to 1P illumination intensity
intVal = intVal./max(intVal(:));

refAsk = (slmXYZ(1:3, 1:end-holdback))';
refGet = intVal(1:end-holdback);

SLMtoPower =  polyfitn(refAsk,refGet,modelterms);

Ask = (slmXYZ(1:3, end-holdback:end))';
True = intVal(end-holdback:end)';

Get = polyvaln(SLMtoPower,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);

figure(1);clf
subplot(2,3,1)
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],intVal,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Measured Power (converted to 1p)')
colorbar
axis square

subplot(2,3,2)
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)'),'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated Power Norm.')
colorbar
axis square

subplot(2,3,4)
% scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal','filled');
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],basVal,'filled');

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Raw Fluorescence')
colorbar
axis square

subplot(2,3,5)
c = sqrt((polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal').^2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error RMS (A.U.)')
colorbar
axis square

subplot(2,3,3)
c = (polyvaln(SLMtoPower,slmXYZ(1:3,:)').^2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated 2P Power')
colorbar
axis square

subplot(2,3,6)
normVal = basVal./max(basVal(:));

c = (polyvaln(SLMtoPower,slmXYZ(1:3,:)').^2)-normVal';
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error 2P Power')
colorbar
axis square


disp(['The RMS error: ' num2str(meanRMS) ' A.U. Power Estimate']);
disp(['The Max power error: ' num2str(max(RMS)*100) '% of request']);

CoC.SLMtoPower = SLMtoPower;
out.CoC = CoC;
out.powerFitmodelTerms = modelterms;

%% THIS BEGINS THE NEW SECTION 3/15/21

% we will use the model to generate a lot of points that we have a pretty
% good guess about where they are located, so we can optimize and multiplex
% the hologram search, using these CoC lets create holograms that shoot a
% pattern into a field of view

%% Determine ScanImage FOV boundaries
% this will also creates the burn pattern (but does not compute the CoC
% for the point yets)

% params
SImatchThreshold = 0; % threshold for being in SI FOV (set to 0 to take whole range)
fineBufferMargin = 0.1;
burnBufferMargin = 0.05;
nBurnGrid = 8; %number of points in the burn grid
FractionOmit = 0.1; % drop points from burn grid


binarySI = ~isnan(SIpeakVal);
SImatchProb = mean(binarySI'); % probability that point was detected aka that was in SI range
SImatchXY = camXYZ(1:2,1:625); % location of points in XYZ

figure(8)
clf
scatter(SImatchXY(1,:),SImatchXY(2,:),[],SImatchProb,'filled');

SIx = SImatchXY(1,SImatchProb>SImatchThreshold);
SIy = SImatchXY(2,SImatchProb>SImatchThreshold);

SIboundary = boundary(SIx',SIy');
SIxboundary = SIx(SIboundary);
SIyboundary = SIy(SIboundary);
hold on
plot(SIxboundary,SIyboundary)

sz = size(bgd);

% for fine targets
SImatchRangeXforFine = [max(min(SIx(SIboundary))-sz(1)*fineBufferMargin,1) ...
    min(max(SIx(SIboundary))+sz(1)*fineBufferMargin,sz(1))];

SImatchRangeYforFine = [max(min(SIy(SIboundary))-sz(2)*fineBufferMargin,1) ...
    min(max(SIy(SIboundary))+sz(2)*fineBufferMargin,sz(2))];
r = rectangle('position',...
    [SImatchRangeXforFine(1) SImatchRangeYforFine(1) SImatchRangeXforFine(2)-SImatchRangeXforFine(1) SImatchRangeYforFine(2)-SImatchRangeYforFine(1)]);
r.EdgeColor='g';


% for hole burning
SImatchRangeX = [max(min(SIx(SIboundary))-sz(1)*burnBufferMargin,1) ...
    min(max(SIx(SIboundary))+sz(1)*burnBufferMargin,sz(1))];

SImatchRangeY = [max(min(SIy(SIboundary))-sz(2)*burnBufferMargin,1) ...
    min(max(SIy(SIboundary))+sz(2)*burnBufferMargin,sz(2))];

r = rectangle('position',...
    [SImatchRangeX(1) SImatchRangeY(1) SImatchRangeX(2)-SImatchRangeX(1) SImatchRangeY(2)-SImatchRangeY(1)]);
r.EdgeColor='r';

rline = line(NaN,NaN,'LineWidth',1','LineStyle', '-','color','r');
gline = line(NaN,NaN,'LineWidth',1','LineStyle', '-','color','g');

colorbar
legend('Prob of SI FOV','Detected SI FOV','Burn Boundary Box','Calib Boundary Box')
xlabel('Camera X pixels')
ylabel('Camera Y Pixels')
% axis equal
axis image

xpts = linspace(SImatchRangeX(1),SImatchRangeX(2),nBurnGrid);
ypts = linspace(SImatchRangeY(1),SImatchRangeY(2),nBurnGrid);

XYpts = [];
for i=1:nBurnGrid
    for k=1:nBurnGrid
        XYpts(:,end+1) = [xpts(i) ypts(k)];
    end
end

XYptsInBounds = inpolygon(XYpts(1,:),XYpts(2,:),SIxboundary,SIyboundary);
XYpts = XYpts(:,XYptsInBounds);

% figure out spacing of pts across zs
zsToBlast = zsToUse; % match to SI Calib
interXdist = xpts(2)-xpts(1);
interYdist = ypts(2)-ypts(1);

gridSide = ceil(sqrt(numel(zsToBlast)));
xOff = round(interXdist/gridSide);
yOff = round(interYdist/gridSide);

% Turn into a more unique looking pattern
numPts = size(XYpts,2);
XYpts(:,randperm(numPts,round(numPts*FractionOmit)))=[];
XYpts = reshape(XYpts,[2 numel(XYpts)/2]);

disp([num2str(size(XYpts,2)) ' points per plane selected. ' num2str(size(XYpts,2)*numel(zsToBlast)) ' total'])

%% Simulate and create new Fine Points
% do a CoC to get more points to shoot

denseFineTimer = tic;

nSimulatedTargs = 10000;
multiholosize = 20;
planes = 7;
holosperplane = 10;

% get basler targets to shoot in from basler range
% for X
a = min(basXYZBackup(1,:));
b = max(basXYZBackup(1,:));
% a = min(basXYZ(1,:));
% b = max(basXYZ(1,:));
% a = 100;
% b = 1000;
% a =SImatchRangeXforFine(1);
% b = SImatchRangeXforFine(2);
r = (b-a).*rand(nSimulatedTargs,1) + a;
rX = round(r);

% for Y
a = min(basXYZBackup(2,:));
b = max(basXYZBackup(2,:));
% a = min(basXYZ(2,:));
% b = max(basXYZ(2,:));
% a = 100;
% b = 1000;
r = (b-a).*rand(nSimulatedTargs,1) + a;
rY = round(r);

% for Z
a = min(basXYZBackup(3,:));
b = max(basXYZBackup(3,:));
% a = min(basXYZ(3,:));
% b = max(basXYZ(3,:));
% a = 5;
% b = 80;
r = (b-a).*rand(nSimulatedTargs,1) + a;
rZ = round(r);

bas2shoot = [rX rY rZ];

testSLM = function_Eval3DCoC(camToSLM, bas2shoot);
expectBas = function_Eval3DCoC(SLMtoCam, testSLM);

testSLM(:,4) = ones(size(testSLM,1),1);

% make sure the SLM vals are within range
excludeMe = testSLM(:,1) < 0 | testSLM(:,1) > 1;
excludeMe = excludeMe | testSLM(:,2) < 0 | testSLM(:,2) > 1;

testSLM = testSLM(~excludeMe,:);
expectBas = expectBas(~excludeMe,:);

% generate multi-target holos that are spread apart
ntotalPoints = multiholosize * planes * holosperplane;
disp(['Using ' num2str(ntotalPoints) ' points in round 2.'])

slm_coords = {};
bas_coords = {};

% group
[idx, ~] = kmeans(expectBas(:,3), planes);

for i=1:planes
    iter=0;
    h=0;
    while 1
        while h < holosperplane
            iter = iter+1;
            if iter > 10000
                disp(['****BAD WARNING! Exited hologram determination loop early. Could not find a suitable hologram for plane ' num2str(i) '.****'])
                break
            end

            targs_this_plane = find(idx==i);
            % choose rand holos
            holo_idxs = randperm(length(targs_this_plane),multiholosize);
            dist = pdist2(expectBas(holo_idxs,:),expectBas(holo_idxs,:));
            temp = rand(size(dist,1));
            dist(find(diag(diag(temp))))=nan;
            if any(dist<100)
                continue
            end
            h = h+1;
            bas_coords{i}{h} = expectBas(targs_this_plane(holo_idxs),:);
            slm_coords{i}{h} = testSLM(targs_this_plane(holo_idxs),:);

            idx(targs_this_plane(holo_idxs)) = -i; %prevent shooting the same target twice. if there are too few in the simulation will error
        end
        break
    end
end



figure(1579)
clf
cmap = colormap(viridis(numel(slm_coords)*holosperplane));
c = 0;
for i = 1:numel(slm_coords)
    hold on
    for j = 1:numel(slm_coords{i})
        c = c + 1;
        hold on
        subplot(1,2,1)
        scatter3(bas_coords{i}{j}(:,1),bas_coords{i}{j}(:,2),bas_coords{i}{j}(:,3), [], cmap(c,:), 'filled')%, 'MarkerFaceAlpha',0.7)
        hold on
        title('Bas Coords')
        subplot(1,2,2)
        scatter3(slm_coords{i}{j}(:,1),slm_coords{i}{j}(:,2),slm_coords{i}{j}(:,3), [], cmap(c,:), 'filled')%, 'MarkerFaceAlpha',0.7)
        title('SLM Coords')
    end
end
%% compute holos

Setup.useGPU = 1;
c = 0;
for i=1:length(slm_coords)
    for j=1:length(slm_coords{i})
        c = c + 1;
        ht = tic;
        disp(['Compiling multi-target hologram ' num2str(c)])

        thisCoord = slm_coords{i}{j};
        thisCoord(:,3) = round(thisCoord(:,3),4); % Added 3/15/21 by Ian for faster compute times

        [ mtholo, ~, ~ ] = function_Make_3D_SHOT_Holos(Setup,thisCoord);
        holos2shoot{i}{j} = mtholo;
        disp(['done in ' num2str(toc(ht)) 's.'])
    end
end

disp('now to shooting...')

%% now repeat multi-target search with new holograms
clear peakValue4 peakDepth4 peakFWHM4 dataUZ4

%% capture holograms

box_range = 20; % distance threshold is set to 100
nframesCapture = 5;

disp('shooting!')

planes = numel(slm_coords);
for i = 1:planes
    holos_this_plane = numel(slm_coords{i});
    disp(['Plane ' num2str(i) ' of ' num2str(planes)])

    % find the mean z of the holo targets and set a range around it
    % rearanging so it doesn't move unnescessarily
    meanz = mean(cellfun(@(x) mean(x(:,3)),bas_coords{i}));
    fineUZ = linspace(meanz-fineRange, meanz+fineRange, finePts);
    dataUZPlane = nan([size(bgd(:,:,1)) finePts holos_this_plane]);

    % for every sutter z plane
    fprintf('Depth: ')
    figure(4);clf;

    for k = 1:finePts
        fprintf([num2str(round(fineUZ(k))) ' ']);

        % move the sutter
        sutter.moveZ(fineUZ(k))

        if i==1
            pause(mvLongWait)
        else
            pause(mvShortWait);
        end

        a = floor(sqrt(holos_this_plane));
        b = ceil(holos_this_plane/a);
        ro = min([a b]);
        co = max([a b]);

        for j= 1:holos_this_plane
            multi_pwr = size(slm_coords{i}{j},1) * pwr;

            slm.feed(holos2shoot{i}{j})

            laser.set_power(multi_pwr)
            data = bas.grab(nframesCapture);
            laser.set_power(0)

            % postprocess
            frame = mean(data,3);
            frame =  max(frame-bgd,0);
            frame = imgaussfilt(frame,2);

            % store into dataUZ(x,y,z-plane)
            dataUZPlane(:,:,k,j) =  frame;

            figure(4)
            subplot(ro,co,j)
            imagesc(frame)
            colorbar
            caxis([0 15]);
            title({['Live Data. Depth ' num2str(round(fineUZ(k)))] ; ['Plane: ' num2str(i) '. Set ' num2str(j)]})
            drawnow

        end
    end
    fprintf('n')

    for j= 1:holos_this_plane
        fineUZ4{i}{j} = fineUZ;
        dataUZ4{i}{j} = uint8(dataUZPlane(:,:,:,j)); %brought out of for loop

        dataUZ = dataUZPlane(:,:,:,j);

        % move sutter back to reference
        sutter.moveToRef()
        pause(mvShortWait)

        % target parsing, might do later instead
        for targ = 1:size(slm_coords{i}{j},1)

            expected_xyz = bas_coords{i}{j}(targ,:);
            [x, y] = size(bgd);

            targX = round(expected_xyz(1)-box_range:expected_xyz(1)+box_range);
            targY = round(expected_xyz(2)-box_range:expected_xyz(2)+box_range);

            if max(targX)>x
                targX = expected_xyz(1)-box_range:x;
            end
            if max(targY)>y
                targY = expected_xyz(2)-box_range:y;
            end
            if min(targX)<1
                targX = 1:expected_xyz(1)+box_range;
            end
            if min(targY)<1
                targY = 1:expected_xyz(2)+box_range;
            end

            try
                % method 1 - rely on XY from first step
                targ_stack = squeeze(max(max(dataUZ(targX,targY,:))));
                mxProj = max(dataUZ(targX,targY,:),[],3);
                [ holo_x,holo_y ] =function_findcenter(mxProj );
                xyFine4{i}{j}(:,targ) = [holo_x+(min(targX)), holo_y+(min(targY))];
            catch
                targ_stack = nan(finePts,1);
                xyFine4{i}{j}(:,targ) =[nan, nan];
            end


            try
                ff = fit(fineUZ', targ_stack, 'gauss1');
                peakValue4{i}{j}(targ) = ff.a1;
                peakDepth4{i}{j}(targ) = ff.b1;
                peakFWHM4{i}{j}(targ) = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
            catch
                disp(['Error on fit! Holo: ', num2str(j), ' Target: ', num2str(targ)])
                peakValue4{i}{j}(targ) = NaN;
                peakDepth4{i}{j}(targ) = NaN;
                peakFWHM4{i}{j}(targ) = NaN;
            end
        end
    end
end

fineT = toc(denseFineTimer);
disp(['Dense Fine Fits took ' num2str(fineT) 's']);

%% New Fits with denser Fine
denseFitsTimer = tic;

c=0;
slmXYZextra = [];
baxXYZextra =[];
basValextra=[];
FWHMValExtra = [];
peakDepthValExtra =[];

for i=1:planes
    for j=1:holos_this_plane
        for targ = 1:size(slm_coords{i}{j},1)
            c=c+1;
            slmXYZextra(c,:) = slm_coords{i}{j}(targ,:);
            baxXYZextra(c,:) = xyFine4{i}{j}(:,targ);
            basValextra(c) = peakValue4{i}{j}(targ);
            FWHMValExtra(c) = peakFWHM4{i}{j}(targ);
            peakDepthValExtra(c) = peakDepth4{i}{j}(targ);
        end
    end
end

%% exclude trials

slmXYZ4 = slmXYZextra';
basXYZ4 = [baxXYZextra peakDepthValExtra']';
basVal4 = basValextra;
FWHMVal4 = FWHMValExtra;%added 7/20/2020 -Ian

excludeTrials = all(basXYZ4(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

% excludeTrials = excludeTrials | basVal4>260; %max of this camera is 255

basDimensions = size(bgd);
excludeTrials = excludeTrials | basXYZ4(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ4(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ4(3,:)<-25; %9/19/19 Ian Added to remove systematic low fits
excludeTrials = excludeTrials | basXYZ4(3,:)>150;


excludeTrials = excludeTrials | any(isnan(basXYZ4(:,:)));
excludeTrials = excludeTrials | basVal4<1; 
excludeTrials = excludeTrials | basVal4>6;

slmXYZBackup = slmXYZ4(:,~excludeTrials);
basXYZBackup = basXYZ4(:,~excludeTrials);
basValBackup = basVal4(:,~excludeTrials);
FWHMValBackup = FWHMVal4(~excludeTrials);

disp(['Number of trials excluded: ' num2str(sum(excludeTrials))])

figure(1922)
clf
subplot(1,2,1)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 75, basValBackup, 'Filled')
colorbar
colormap default
title({'Second Denser Fine';'basXYZ and basVals (fine)'})

subplot(1,2,2)
hold on
basx = 1:size(basVal4,2);
plot(basx, basVal4,'ko')
plot(basx(~excludeTrials), basVal4(~excludeTrials), 'o')
legend('Excluded','Included')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')

%% fit SLM to Camera
%use model terms

errScalar = 2.5;
holdback = 500;
pxPerMu = 1;


basXYZ4 = basXYZBackup;
slmXYZ4 = slmXYZBackup;
basVal4 = basValBackup;

disp('Fitting SLM to Camera')
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1 ; 1 1 1 ;...
    2 0 0; 0 2 0; 0 0 2;  ...
    2 0 1; 2 1 0; 0 2 1; 1 2 0; 0 1 2;  1 0 2; ... ];  %XY spatial calibration model for Power interpolations
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2;];

reOrder = randperm(size(slmXYZ4,2));
slmXYZ4 = slmXYZ4(:,reOrder);
basXYZ4 = basXYZ4(:,reOrder);


refAsk = (slmXYZ4(1:3,1:end-holdback))';
refGet = (basXYZ4(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

figure(1977)
clf
subplot(1,2,1)
ax = gca();
[SLMtoCam, ~] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,ax);
title('SLM to Cam v2')

Ask = refAsk;
True = refGet;
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(103);
clf
subplot(1,2,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
% legend('Measured targets', 'Estimated Targets');
title({'Reference Data'; 'SLM to Camera'})

refRMS = sqrt(sum((Get-True).^2,2));
subplot(1,2,2)
scatter3(True(:,1),True(:,2),True(:,3),[],refRMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title({'Reference Data'; 'RMS Error in position'})
caxis([0 30])


Ask = (slmXYZ4(1:3,end-holdback:end))';
True = (basXYZ4(1:3,end-holdback:end))';
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(101);
clf
subplot(1,3,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')


ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
legend('Measured targets', 'Estimated Targets');
title('SLM to Camera')

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp('Error based on Holdback Data...')
disp(['The RMS error: ' num2str(meanRMS) ' pixels for SLM to Camera']);


disp(['Thats approx ' num2str(meanRMS/pxPerMu) ' um']);

xErr = sqrt(sum((Get(:,1)-True(:,1)).^2,2));
yErr = sqrt(sum((Get(:,2)-True(:,2)).^2,2));
zErr = sqrt(sum((Get(:,3)-True(:,3)).^2,2));

disp('Mean:')
disp(['X: ' num2str(mean(xErr)/pxPerMu) 'um. Y: ' num2str(mean(yErr)/pxPerMu) 'um. Z: ' num2str(mean(zErr)) 'um.']);
disp('Max:')
disp(['X: ' num2str(max(xErr)/pxPerMu) 'um. Y: ' num2str(max(yErr)/pxPerMu) 'um. Z: ' num2str(max(zErr)) 'um.']);

subplot(1,3,2)
scatter3(True(:,1),True(:,2),True(:,3),[],RMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title('RMS Error in position')


refAsk = (basXYZ4(1:3,1:end-holdback))';
refGet = (slmXYZ4(1:3,1:end-holdback))';

%  camToSLM = function_3DCoC(refAsk,refGet,modelterms);
figure(1977)
subplot(1,2,2)
ax = gca();
[camToSLM, ~] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0, ax);
title('Cam to SLM v2')

Ask = (basXYZ4(1:3,end-holdback:end))';
True = (slmXYZ4(1:3,end-holdback:end))';
Get = function_Eval3DCoC(camToSLM,Ask);

figure(101)
subplot(1,3,3)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Depth units')
legend('Measured targets', 'Estimated Targets');
title('Camera to SLM')

CoC.camToSLM=camToSLM;
CoC.SLMtoCam = SLMtoCam;

out.CoC=CoC;
out.CoCmodelterms = modelterms;

rtXYZ = function_Eval3DCoC(SLMtoCam,function_Eval3DCoC(camToSLM,basXYZ4(1:3,end-holdback:end)'));

err = sqrt(sum((rtXYZ - basXYZ4(1:3,end-holdback:end)').^2,2));
meanRTerr = nanmean(err);
disp(['The Mean Round Trip RMS error: ' num2str(meanRTerr) ' pixels (' num2str(meanRTerr/pxPerMu) ' um) camera to SLM to camera']);

%% fit power as a function of SLM
disp('Fitting Power as a function of SLM')

slmXYZ4 = slmXYZBackup;
basVal4 = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];

intVal = basVal4;
intVal = sqrt(intVal); %convert fluorescence intensity (2P) to 1P illumination intensity
intVal=intVal./max(intVal(:));

refAsk = (slmXYZ4(1:3,1:end-holdback))';
refGet = intVal(1:end-holdback);

SLMtoPower =  polyfitn(refAsk,refGet,modelterms);

Ask = (slmXYZ4(1:3,end-holdback:end))';
True = intVal(end-holdback:end)';

Get = polyvaln(SLMtoPower,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);

figure(1111);clf
subplot(2,3,1)
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],intVal,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Measured Power (converted to 1p)')
colorbar
% axis square

subplot(2,3,2)
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],polyvaln(SLMtoPower,slmXYZ4(1:3,:)'),'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated Power Norm.')
colorbar
% axis square

subplot(2,3,4)
% scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal','filled');
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],basVal4,'filled');

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Raw Fluorescence')
colorbar
% axis square

subplot(2,3,5)
c = sqrt((polyvaln(SLMtoPower,slmXYZ4(1:3,:)')-intVal').^2);
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error RMS (A.U.)')
colorbar
% axis square

subplot(2,3,3)
c = (polyvaln(SLMtoPower,slmXYZ4(1:3,:)').^2);
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated 2P Power')
colorbar
% axis square

subplot(2,3,6)
normVal = basVal4./max(basVal4(:));

c = (polyvaln(SLMtoPower,slmXYZ4(1:3,:)').^2)-normVal';
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error 2P Power')
colorbar
% axis square


disp(['The RMS error: ' num2str(meanRMS) ' A.U. Power Estimate']);
disp(['The Max power error: ' num2str(max(RMS)*100) '% of request']);

CoC.SLMtoPower = SLMtoPower;
out.CoC = CoC;
out.powerFitmodelTerms = modelterms;

%% Plot FWHM
FWHM = FWHMValBackup;
depth = basXYZBackup(3,:);
slmXYZ = slmXYZBackup;

figure(1001); clf
subplot(1,2,1);
plot(FWHM,depth,'o')
% plot(FWHM,slmCoords(3,:),'o')

ylabel('Axial Depth \mum')
xlabel('FWHM \mum')
ylim([-25 125])
xlim([7.5 50])

refline(0,0)
refline(0,60)


subplot(1,2,2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],FWHM,'filled')
% scatter3(slmXYZ(1,:),slmXYZ(2,:),depth,[],FWHM,'filled')
caxis([10 50])
h= colorbar;
xlabel('SLM X')
ylabel('SLM Y')
zlabel('SLM Z')
set(get(h,'label'),'string','FWHM \mum')

fprintf(['FWHM in the typical useable volume (0 to 100um) is: ' num2str(mean(FWHM(depth>0 & depth<60))) 'um\n'])


finalFitsT = toc(denseFitsTimer);
%% Plot Hole Burn Stuff
tCompileBurn = tic;

%figure(6);
%scatter(XYpts(1,:),XYpts(2,:),'o');

figure(44); clf

clear XYtarg SLMtarg
for i = 1:numel(zsToBlast)
    a = mod(i-1,gridSide);
    b = floor((i-1)/gridSide);

    XYuse = bsxfun(@plus,XYpts,([xOff*a yOff*b])');
    optoZ = zsToBlast(i);

    zOptoPlane = ones([1 size(XYuse,2)])*optoZ;

    Ask = [XYuse; zOptoPlane];
    estCamZ = polyvaln(OptZToCam,Ask');
    meanCamZ(i) = nanmean(estCamZ); %for use by sutter
    Ask = [XYuse; estCamZ'];
    estSLM = function_Eval3DCoC(camToSLM,Ask');
    estPower = polyvaln(SLMtoPower,estSLM);

    % negative DE restrictions
    ExcludeBurns = ((estPower<=0) | (estSLM(:,1)<slmXrange(1)) | (estSLM(:,1)>slmXrange(2)) | (estSLM(:,2)<slmYrange(1)) | (estSLM(:,2)>slmYrange(2))); %don't shoot if you don't have the power
    estSLM(ExcludeBurns,:)=[];
    estPower(ExcludeBurns)=[];
    XYuse(:,ExcludeBurns)=[];
    zOptoPlane(ExcludeBurns)=[];
    estCamZ(ExcludeBurns)=[];

    XYtarg{i} = [XYuse; zOptoPlane];
    SLMtarg{i} = [estSLM estPower];

    subplot(1,2,1)
    scatter3(XYuse(1,:),XYuse(2,:),estCamZ,[],estPower,'filled')

    hold on
    subplot (1,2,2)
    scatter3(estSLM(:,1),estSLM(:,2),estSLM(:,3),[],estPower,'filled')

    disp([num2str(min(estSLM(:,3))) ' ' num2str(max(estSLM(:,3)))])
    hold on
end

subplot(1,2,1)
title('Targets in Camera Space')
zlabel('Depth \mum')
xlabel('X pixels')
ylabel('Y pixels')

subplot(1,2,2)
title('Targets in SLM space')
xlabel('X SLM')
ylabel('Y SLM')
zlabel('Z SLM')
c = colorbar;
c.Label.String = 'Estimated Power';

%% Burn Holes
disp('Compiling Holos To Burn')

blankHolo = zeros(1024, 1024);

clear tempHololist
for k = 1:numel(zsToBlast)
    parfor i=1:size(XYtarg{k},2)
        t=tic;
        fprintf(['Compiling Holo ' num2str(i) ' for depth ' num2str(k)]);
        subcoordinates =  [SLMtarg{k}(i,1:3) 1];
        %check to avoid out of range holos; added 11/1/19
        %12/5/19 now it allows negative zs
        if ~any(subcoordinates(1:2)>1 | subcoordinates(1:2) <0)
            DE(i) = SLMtarg{k}(i,4);
            [ Hologram,~,~ ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates );
        else
            DE(i)= 0;
            Hologram = blankHolo;
        end
        tempHololist(:,:,i)=Hologram;
        fprintf([' took ' num2str(toc(t)) 's\n']);
    end
    holos{k}=tempHololist;
    Diffraction{k}=DE;
end
disp(['Compiling Done took ' num2str(toc(tCompileBurn)) 's']);

compileBurnT = toc(tCompileBurn);

%% do the burning

% params
si_power = 12; % percent
burnPowerMultiplier = 8;
burnTime = 0.5; % in seconds, very rough and not precise
blastPowerCap = 2; % Watts

disp('Blasting Holes for SI to SLM alignment, this will take about an hour and take 25Gb of space')
tBurn = tic;

%confirm that SI computer in eval Mode
mssend(SISocket,'1+2');
invar=[];
while ~strcmp(num2str(invar),'3') %stupid cludge so that [] read as false
    invar = msrecv(SISocket,0.01);
end
disp('linked')

%setup acquisition

numVol = 7; %number of SI volumes to average
baseName = '''calib''';

mssend(SISocket,['hSI.hStackManager.arbitraryZs = [' num2str(zsToBlast) '];']);
mssend(SISocket,['hSI.hStackManager.numVolumes = [' num2str(numVol) '];']);
mssend(SISocket,'hSI.hStackManager.enable = 1 ;');
mssend(SISocket,'hSI.hBeams.pzAdjust = 0;');
mssend(SISocket,['hSI.hBeams.powers = ' num2str(si_power) ';']);
mssend(SISocket,'hSI.extTrigEnable = 0;'); %savign
mssend(SISocket,'hSI.hChannels.loggingEnable = 1;'); %savign
mssend(SISocket,'hSI.hScan2D.logFilePath = ''D:\Calib\Temp'';');
mssend(SISocket,['hSI.hScan2D.logFileStem = ' baseName ';']);
mssend(SISocket,'hSI.hScan2D.logFileCounter = 1;');
mssend(SISocket,'hSICtl.updateView;');

%clear invar
invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = msrecv(SISocket,0.01);
end

mssend(SISocket,'30+7');
invar=[];
while ~strcmp(num2str(invar),'37')
    invar = msrecv(SISocket,0.01);
end
disp('completed parameter set')

%%Burn

%AcquireBaseline
disp('Acquire Baseline')

mssend(SISocket,'hSI.startGrab()');
invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = msrecv(SISocket,0.01);
end
wait = 1;
while wait
    mssend(SISocket,'hSI.acqState;');
    invar = msrecv(SISocket,0.01);
    while isempty(invar)
        invar = msrecv(SISocket,0.01);
    end

    if strcmp(invar,'idle')
        wait=0;
        disp('Ready for Next')
    else
        %             disp(invar)
    end
end



disp('Now Burning')

for k=1:numel(zsToBlast)

    offset = round(meanCamZ(k));
    sutter.moveZ(offset)

    if k==1
        pause(mvLongWait)
    else
        pause(mvShortWait);
    end

    tempHololist=holos{k};

    for i=1:size(XYtarg{k},2)%1:size(XYuse,2)
        t=tic;
        fprintf(['Blasting Hole ' num2str(i) '. Depth ' num2str(zsToBlast(k))]);
        slm.feed(tempHololist(:,:,i))

        DE = Diffraction{k}(i);
        % if Diffraction efficiency too low just don't even burn
        if DE<.05 
            DE=inf;
        end
        blastPower = pwr*burnPowerMultiplier/DE;
        
        % blastPower is in mW and and blastPowerCap is in W, so be sure to
        % convert
        if blastPower > blastPowerCap*1000
            blastPower = blastPowerCap*1000;
        end

        stimT=tic;
        laser.set_power(blastPower)
        while toc(stimT)<burnTime
        end
        laser.set_power(0)

        %flush and handshake
        laser.flush()
        %re send 0
        laser.set_power(0)

        mssend(SISocket,'hSI.startGrab()');
        invar = msrecv(SISocket,0.01);
        while ~isempty(invar)
            invar = msrecv(SISocket,0.01);
        end

        wait = 1;
        while wait
            mssend(SISocket,'hSI.acqState');
            invar = msrecv(SISocket,0.01);
            while isempty(invar)
                invar = msrecv(SISocket,0.01);
            end

            if strcmp(invar,'idle')
                wait=0;
                %             disp(['Ready for Next'])
            else
                %             disp(invar)
            end
        end
        disp([' Took ' num2str(toc(t)) 's'])
    end
end

sutter.moveToRef()

burnT = toc(tBurn);
disp(['Done Burning. Took ' num2str(burnT) 's']);

disp('Done with Lasers and ScanImage now, you can turn it off')

%% Move file to modulation

disp('Moving files')
tMov = tic;

%on ScanImage Computer
% destination = '''F:\frankenshare\FrankenscopeCalib''' ;
destination = '''F:\Calib\Temp''';
source = '''D:\Calib\Temp\calib*''';

%clear invar
invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = msrecv(SISocket,0.01);
end


mssend(SISocket,['movefile(' source ',' destination ')']);
invar = msrecv(SISocket,0.01);
while isempty(invar)
    invar = msrecv(SISocket,0.01);
end
disp(['Moved. Took ' num2str(toc(tMov)) 's']);
MovT= toc(tMov);


%% read/compute frame

mssend(SISocket,'end');
%%

tLoad = tic;
% pth = 'W:\frankenshare\FrankenscopeCalib'; %On this compute
% pth ='F:\Temp'; %temp fix b/c frankenshare down
pth = 'F:\Calib\Temp';
files = dir(pth);

baseN = eval(baseName);

[dummy fr] = bigread3(fullfile(pth,files(3).name) );

nOpto = numel(zsToBlast);
nBurnHoles = size(XYtarg{1},2);

baseFr = mean(fr(:,:,1:nOpto:end),3);%mean(fr(:,:,1:nOpto:end),3);%Probably more accurate to just do correct zoom, but sometimes having difficulty

k=1;c=0; SIXYZ =[];
for i=4:numel(files)
    t = tic;
    fprintf(['Loading/Processing Frame ' num2str(i)]);
    try
        [dummy fr] = bigread3(fullfile(pth,files(i).name) );
        if c>=nBurnHoles
            k=k+1;
            c=0;
            nBurnHoles = size(XYtarg{k},2);
        end
        c=c+1;

        Frame = mean(fr(:,:,:),3);%mean(fr(:,:,k:nOpto:end),3); %Probably more accurate to just do correct zoom, but sometimes having difficulty
        Frames{k}(:,:,c) = Frame;

        if c>1
            baseFrame = Frames{k}(:,:,c-1);

            %try to exclude those very bright spots
            maskFR = imgaussfilt(Frame,3) - imgaussfilt(Frame,15);
            mask = maskFR > mean(maskFR(:))+6*std(maskFR(:));

            %remove the low frequency slide illumination differences
            filtNum = 2;
            frameFilt = imgaussfilt(Frame,filtNum);
            baseFilt = imgaussfilt(baseFrame,filtNum);


            toCalc = (baseFrame-baseFilt) - (Frame-frameFilt);
            toCalc(mask)=0;

            %             testFr = Frames{k}(:,:,c-1) - Frame;
            [ x,y ] =function_findcenter(toCalc);

            figure(333)
            clf
            subplot(1,3,1)
            imagesc(Frame)

            subplot(1,3,2)
            imagesc(frameFilt)

            subplot(1,3,3)
            imagesc(toCalc)
            hold on
            scatter(y,x,[],'r')
            %             pause
        else
            x = 0;
            y=0;
        end
    catch
        fprintf('\nError in Hole analysis... probably loading.')
        x = 0;
        y=0;
    end


    SIXYZ(:,end+1) = [x,y,zsToBlast(k)];
    disp([' Took ' num2str(toc(t)) ' s']);
end

SIXYZbackup=SIXYZ;
disp(['Done Loading/Processing SI files. Took ' num2str(toc(tLoad)) 's'])
loadT = toc(tLoad);


%% do non-cv SI to cam calculation
errScalar = 2.5;

burnFitsTimer = tic;

modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

cam3XYZ = [XYtarg{:};];
SIXYZ = SIXYZbackup;

cam3XYZ=cam3XYZ(:,1:size(SIXYZ,2));

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;
cam3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];

refAsk = SIXYZ(1:3,:)';
refGet = (cam3XYZ(1:3,:))';

figure(2594);clf
subplot(1,2,1)
ax = gca();
[SItoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0, ax);
title('SI to Cam')

subplot(1,2,2)
ax = gca();
[CamToSI, trialN] = function_3DCoCIterative(refGet,refAsk,modelterms,errScalar,0, ax);
title('Cam to SI')

CoC.CamToSI = CamToSI;
CoC.SItoCam = SItoCam;
out.CoC=CoC;

%% alternate calculation
errScalar = 2.5;

tempSLM = cellfun(@(x) x',SLMtarg,'UniformOutput',false);
slm3XYZ = [tempSLM{:}];
SIXYZ = SIXYZbackup;

slm3XYZ=slm3XYZ(1:3,1:size(SIXYZ,2));

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;

slm3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];

refAsk = SIXYZ(1:3,:)';
refGet = (slm3XYZ(1:3,:))';

figure(2616);clf
subplot(1,2,1)
ax = gca();
[SItoSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0, ax);
title('SI to SLM')
subplot(1,2,2)
ax = gca();
[SLMtoSI, trialN] = function_3DCoCIterative(refGet,refAsk,modelterms,errScalar,0,ax);
title('SLM to SI')

CoC.SItoSLM = SItoSLM;
CoC.SLMtoSI = SLMtoSI;


%% Calculate round trip errors
numTest = 10000;

rangeX = [0 511];%[0 511];
rangeY = [0 511];%[0 511];
rangeZ = [0 55];% Make Sure to match this to the correct range for this optotune;

clear test;
valX = round((rangeX(2)-rangeX(1)).*rand(numTest,1)+rangeX(1));
valY = round((rangeY(2)-rangeY(1)).*rand(numTest,1)+rangeY(1));
valZ = round((rangeZ(2)-rangeZ(1)).*rand(numTest,1)+rangeZ(1));

test = [valX valY valZ];
%%display
test2 = function_SLMtoSI(function_SItoSLM(test,CoC),CoC);
ER1xy = test2(:,1:2)-test(:,1:2);
RMSE1xy = sqrt(sum(ER1xy'.^2));

SIpxPerMu = 512/800;

ER1z = test2(:,3)-test(:,3);
RMSE1z = abs(ER1z);

meanE1rxy = mean(RMSE1xy);
meanE1rz = mean(RMSE1z);

figure(12);clf
subplot(4,2,1)
histogram(RMSE1xy/SIpxPerMu,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')
title({'4 Step CoC'; ['Mean RMS err: ' num2str(meanE1rxy) '\mum']})

subplot(4,2,2)
histogram(RMSE1z,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE1rz) ' optotune Units'])

estSLM = function_Eval3DCoC(CoC.SItoSLM,test);
test2 = function_Eval3DCoC(CoC.SLMtoSI,estSLM);
ER2xy = test2(:,1:2)-test(:,1:2);
RMSE2xy = sqrt(sum(ER2xy'.^2));

SIpxPerMu = 512/800;

ER2z = test2(:,3)-test(:,3);
RMSE2z = abs(ER2z);
meanE2rxy = mean(RMSE2xy);
meanE2rz = mean(RMSE2z);

subplot(4,2,3)
histogram(RMSE2xy/SIpxPerMu,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')
title({'1 Step CoC'; ['Mean RMS err: ' num2str(meanE2rxy) '\mum']})

subplot(4,2,4)
histogram(RMSE2z,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE2rz) ' optotune Units'])


estSLM = function_Eval3DCoC(CoC.SItoSLM,test);
estSIasym = function_SLMtoSI(estSLM,CoC);

ERA = estSIasym-test;
RMSErAxy = sqrt(sum(ERA(:,1:2)'.^2));
RMSErAz = abs(ERA(:,3));


subplot(4,2,5)
histogram(RMSErAxy,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')

meanE3rxy = mean(RMSErAxy);
meanE3rz = mean(RMSErAz);
title({'Asymetric CoC; 1S Forward, 4S Reverse'; ['Mean RMS err: ' num2str(meanE3rxy) '\mum']})

subplot(4,2,6)
histogram(RMSErAz,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE3rz) ' optotune Units'])



estSLM2 = function_SItoSLM(test,CoC);
estSLM2 = estSLM2(:,1:3);

estSIasym2 = function_Eval3DCoC(CoC.SLMtoSI,estSLM2);


ERA2 = estSIasym2-test;
RMSErAxy = sqrt(sum(ERA2(:,1:2)'.^2));
RMSErAz = abs(ERA2(:,3));


subplot(4,2,7)%aysmetric reverse; foward with 4 chan
histogram(RMSErAxy,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')

meanE3rxy = mean(RMSErAxy);
meanE3rz = mean(RMSErAz);
title({'Asymetric CoC reverse. 4S Forward, 1S Reverse'; ['Mean RMS err: ' num2str(meanE3rxy) '\mum']})

subplot(4,2,8)
histogram(RMSErAz,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE3rz) ' optotune Units'])

%%Plot scatter
N=10000;


figure(13);clf
subplot(1,2,1)
val=RMSErAxy;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated XY error, both methods')

subplot(1,2,2)
val=RMSErAz;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated Z error, both methods')




figure(600);clf
subplot(1,2,1)
val=RMSE1xy/SIpxPerMu;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated XY error, 1st methods')

subplot(1,2,2)
val=RMSE1z;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated Z error, 1st methods')

burnFitsT = toc(burnFitsTimer);

%% Save Output Function
pathToUse = 'C:\Users\Holography\Desktop\SLM_Management\Calib_Data';
disp('Saving...')
tSave = tic;

save(fullfile(pathToUse,[date '_Calib.mat']),'CoC')
save(fullfile(pathToUse,'ActiveCalib.mat'),'CoC')

pth = 'C:\Users\Holography\Desktop\calibs';
save(fullfile(pth,[date '_Calib.mat']),'CoC')
save(fullfile(pth,'ActiveCalib.mat'),'CoC')
save(fullfile(pth,'CalibWorkspace.mat'), '-v7.3');


% times.saveT = toc(tSave);
% times.burnFitsT = burnFitsT;
% times.loadT = loadT;
% times.MovT = MovT;
% times.burnT = burnT;
% times.compileBurnT = compileBurnT;
% times.finalFitsT = finalFitsT; %Final Fits Camera to SLM
% times.finalFineT = fineT; %Second Dense Fine
% times.intermediateFitsT = intermediateFitsT;
% times.intermediateT = multiT; %first pass fine fit
% times.coarseFitsT = fitsT; %Coarse
% times.coarseT = coarseT; %Coarse Fit
% times.siT = siT;
% times.singleCompileT = singleCompileT;
% times.multiCompileT = multiCompileT;
% times.manualT = tManual; %time spent doing manual setup.

totT = toc(tBegin);
times.totT = totT;
save(fullfile(pathToUse,'CalibWorkspace_v73.mat'), '-v7.3');
% save(fullfile(pathToUse,['CalibWorkspace_will_' date '.mat']))
disp(['Saving took ' num2str(toc(tSave)) 's']);

disp(['All Done, total time from begining was ' num2str(toc(tBegin)) 's. Bye!']);

%% use this to run
%[SLMXYZP] = function_SItoSLM(SIXYZ,CoC);


slm.stop()  
bas.stop()
clear bas slm sutter





