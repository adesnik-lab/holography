%% Pathing and setup

% guards against hitting F5 instead of F9
in = input('Run Calibration? (y/n) ','s');
if ~strcmp(in,'y')
    disp('Did not detect ''y'' so not executing')
    return
end

makePaths()

try function_close_sutter( Sutter ); end
try [Setup.SLM ] = Function_Stop_SLM( Setup.SLM ); end

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

if Setup.useGPU
    disp('Getting gpu...'); %this can sometimes take a while at initialization
    g= gpuDevice;
end

[Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
[ Setup.SLM ] = Function_Start_SLM( Setup.SLM );

Setup.Sutterport ='COM3';
try; function_close_sutter( Sutter ); end
[ Sutter ] = function_Sutter_Start( Setup );

bas = bascam();
bas.start()

% these may no longer be needed...
castImg = bas.castFun;
castAs = bas.castAs;
camMax = bas.camMax;

disp('Ready')
%% look for objective in 1p

bas.preview()

%% connect to DAQ computer

%run this first then code on daq
disp('Waiting for msocket communication From DAQ... ')
%then wait for a handshake
srvsock = mslisten(42116);
masterSocket = msaccept(srvsock,15);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);

invar = [];

while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
end
fprintf('done.\r')
%% connect to ScanImage computer

% run this code first, then 'autoCalibSI' on SI computer
fprintf('Waiting for msocket communication to ScanImage Computer... ')
% then wait for a handshake
srvsock2 = mslisten(42039);
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

pwr = 10;
slmCoords = [.4 .45 -.05 1];

disp(['Individual hologram power set to ' num2str(pwr) 'mW.'])

DEestimate = DEfromSLMCoords(slmCoords);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
Function_Feed_SLM(Setup.SLM, Holo);
mssend(masterSocket, [pwr/1000 1 1]);
bas.preview()
mssend(masterSocket, [0 1 1]);

%% center FOV

disp('Find Focal Plane, Center and Zero the Sutter')
disp('Leave the focus at Zoom 1. at a power that is less likely to bleach (14% 25mW)') %25% 8/16/19
disp('Don''t forget to use Ultrasound Gel on the objective so it doesn''t evaporate')
mssend(masterSocket,[0 1 1]);

bas.preview()

input('Turn off Focus and press any key to continue');

Sutter.Reference = getPosition(Sutter.obj);

mssend(SISocket,[0 0]);

disp('Make Sure the DAQ computer is running testMultiTargetsDAQ. and the SI computer running autoCalibSI');
disp('also make those names better someday')
disp('Make sure both lasers are on and the shutters open')
disp('Scanimage should be idle, nearly in plane with focus. and with the gain set high enough to see most of the FOV without saturating')

position = Sutter.Reference;
position(3) = position(3) + 100;
moveTo(Sutter.obj, position);
disp('testing the sutter double check that it moved to reference +100');
input('Ready to go (Press any key to continue)');

moveTo(Sutter.obj, Sutter.Reference);

tManual = toc(tBegin);

%% Create initial holos for coarse search

npts = 200;

slmXrange = [0.02 0.98];
slmYrange = [0.02 0.98];
slmZrange = [-.07 0.10];

slmCoordsInitial = generateRandomHolos(slmXrange, slmYrange, slmZrange, npts);

figure(1);
scatter3(slmCoordsInitial(1,:), ...
         slmCoordsInitial(2,:), ...
         slmCoordsInitial(3,:))
title('Initial Holograms')

fprintf('Compiling initial holograms... ')
hololist = compileSingleTargetHolos(slmCoordsInitial);
fprintf('done.\r')

%% Measure background noise in camera

bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

meanBgd = mean(bgd_frames, 'all');
stdBgd = std(bgd_frames, [], 'all');

%% Capture ScanImage/Optotune Planes Data

clear SIdepthData

optotunePlanes  = 0:5:75;
sutterPlanes    = -25:5:150;
nframesCapture  = 5;

sz = size(bgd);
gridpts = 25;
xs = round(linspace(1,sz(1),gridpts+2));
ys = round(linspace(1,sz(2),gridpts+2));
xs([1 end])=[];
ys([1 end])=[];

[X,Y] = meshgrid(xs, ys);
XYSI = [X(:) Y(:)]';

clear dimx dimy XYSI2
range =15;
c=0;
for i=1:gridpts
    for k=1:gridpts
        c=c+1;
        dimx(:,c) = xs(i)-range:xs(i)+range;
        dimy(:,c) = ys(k)-range:ys(k)+range;
        
        XYSI2(:,c) = [xs(i) ys(k)];
    end
end

disp(['We will collect ' num2str(numel(zsToUse)) ' planes.'])

SIVals = zeros([SIpts c numel(zsToUse)]);

% for each optotune plane, locate in basler with sutter Z
for k=1:numel(optotunePlanes)
    z = optotunePlanes(k);

    fprintf(['Testing plane ' num2str(z) ': ']);
    
    mssend(SISocket,[z 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(SISocket,0.01);
    end
    
    dataUZ = zeros([sz SIpts]);

    % move the sutter through each optotune plane
    for i = 1:numel(sutterPlanes)
        fprintf([num2str(round(sutterPlanes(i))) ' ']);
        
        position = Sutter.Reference;
        position(3) = position(3) + (SIUZ(i));
        moveTo(Sutter.obj,position);

        if i==1
            pause(1)
        else
            pause(0.1);
        end
        
        data = bas.grab(nframesCapture);

        data = mean(data, 3);
        frame = max(frame-bgd, 0);
        frame = imgaussfilt(frame, 2);

        dataUZ(:,:,i) =  frame;

    end

    SIdepthImages{k} = dataUZ;

    moveTo(Sutter.obj, Sutter.Reference);
    pause(1)
    
    mssend(SISocket,[z 0]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(SISocket,0.01);
    end
    
end

mssend(SISocket,'end');

disp(['Scanimage calibration done whole thing took ' num2str(toc(tSI)) 's']);

%% Extract and do fits

for k = 1:numel(SIdepthImages)
    for i =1:c
        tmp = SIdepthImages{k};
        SIVals(:,i,k) = squeeze(mean(mean(tmp(dimx(:,i),dimy(:,i),:))));
    end
end

SIfits = extractScanImageFits(SIVals, zsToUse, SIUZ);

% exclusions
SIThreshHold = 3*stdBgd/sqrt(nBackgroundFrames + framesToAcquire);
%  SIThreshHold = 0.2;
excl = SIpeakVal<SIThreshHold;

disp([num2str(numel(SIpeakDepth)) ' points total before exclusions'])
disp([num2str(sum(excl(:))) ' points excluded b/c below threshold'])

SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;

excl = SIpeakDepth<-25 | SIpeakDepth>150;

disp([num2str(sum(excl(:))) ' points excluded b/c too deep'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;
disp([num2str(sum(~isnan(SIpeakDepth), 'all')) ' points remaining'])



















