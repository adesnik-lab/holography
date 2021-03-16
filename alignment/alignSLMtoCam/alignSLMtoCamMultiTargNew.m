function alignSLMtoCamMultiTargNew
%% Pathing
in = input('Run Calibration? (y/n)','s');
if ~strcmp(in,'y')
    disp('Did not detect ''y'' so not executing')
    return
end

try function_close_sutter( Sutter ); end
try function_stopBasCam(Setup); end
try [Setup.SLM ] = Function_Stop_SLM( Setup.SLM ); end


clear;close all;clc
%%
tBegin = tic;
addpath(genpath('C:\Users\Holography\Documents\MATLAB\msocket\'));
rmpath(genpath('C:\Users\Holography\Documents\GitHub\SLM-Managment\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\New_SLM_Code\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\NOVOCGH_Code\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\Basler\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\IanTestCode\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\Will test code\'));
disp('done pathing')

%% Setup Stuff
disp('Setting up stuff...');

[Setup ] = function_loadparameters2();
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

try function_stopBasCam(Setup); end
[Setup] = function_startBasCam(Setup);
disp('Ready')

%% look for objective in 1p

function_BasPreview(Setup);


%% Make mSocketConnections with DAQ and SI Computers

disp('Waiting for msocket communication From DAQ')
%then wait for a handshake
srvsock = mslisten(42167);
masterSocket = msaccept(srvsock,15);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);
%MasterIP = '128.32.177.217';
%masterSocket = msconnect(MasterIP,3002);
% ip = '0.0.0.0'; % accept from any ip
% port = 4444;
% masterSocket = tcpip(ip, port, 'NetworkRole', 'server');
% fopen(masterSocket);


invar = [];

while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
end
disp('communication from Master To Holo Established');
%%
disp('Waiting for msocket communication to ScanImage Computer')
%then wait for a handshake
srvsock2 = mslisten(4239);
SISocket = msaccept(srvsock2,15);
msclose(srvsock2);
sendVar = 'A';
mssend(SISocket, sendVar);
%MasterIP = '128.32.177.217';
%masterSocket = msconnect(MasterIP,3002);

invar = [];

while ~strcmp(invar,'B');
    invar = msrecv(SISocket,.1);
end;
disp('communication from Master To SI Established');

%% Put all Manual Steps First so that it can be automated

%% Set Power Levels

pwr = 17; %updated 3/10/21 for 2 MHz % something like this for 100 divided mode 9/22/20 %40; %13 at full; 50 at 15 divided %70 mW at 100 divided mode 10-29-19
disp(['individual hologram power set to ' num2str(pwr) 'mW']);
%%
disp('Find the spot and check if this is the right amount of power')
slmCoords = [.45 .45 -.05 1]%[0.45 0.45 0 1];
DEestimate = DEfromSLMCoords(slmCoords); %
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[ Holo,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );

blankHolo = zeros([1920 1152]);
% Function_Feed_SLM( Setup.SLM, blankHolo);

Function_Feed_SLM( Setup.SLM, Holo);

mssend(masterSocket,[pwr/1000 1 1]);

function_BasPreview(Setup);
mssend(masterSocket,[0 1 1]);

%% Make Sure you're centered
disp('Find Focal Plane, Center and Zero the Sutter')
disp('Leave the focus at Zoom 1. at a power that is less likely to bleach (14% 25mW)') %25% 8/16/19
disp('Don''t forget to use Ultrasound Gel on the objective so it doesn''t evaporate')
mssend(masterSocket,[0 1 1]);

% function_Basler_Preview(Setup, 5);
function_BasPreview(Setup);

temp = input('Turn off Focus and press any key to continue');
Sutter.Reference = getPosition(Sutter.obj);


mssend(SISocket,[0 0]);

disp('Make Sure the DAQ computer is running testMultiTargetsDAQ. and the SI computer running autoCalibSI');
disp('also make those names better someday')
disp('Make sure both lasers are on and the shutters open')
disp('Scanimage should be idle, nearly in plane with focus. and with the gain set high enough to see most of the FOV without saturating')


position = Sutter.Reference;
position(3) = position(3) + 100;
moveTime=moveTo(Sutter.obj,position);
disp('testing the sutter double check that it moved to reference +100');
temp = input('Ready to go (Press any key to continue)');

position = Sutter.Reference;
moveTime=moveTo(Sutter.obj,position);

%% Create a random set of holograms or use flag to reload
disp('First step Acquire Holograms')
reloadHolos =0;
tSingleCompile = tic;
 
if ~reloadHolos
    disp('Generating New Holograms...')
    disp('Everything after this should be automated so sitback and enjoy')
    
    npts = 200; %You can almost get through 750 with water before it evaporates.
    
    RX = 0.4;
    RY = 0.4;
    
    %ranges set by exploration moving holograms looking at z1 fov.
    %updated 3/10/21 (new SLM)
    slmXrange = [0.2 0.85];%9/19/19 [.2 .9]; %[0.125 0.8]; %[0.5-RX 0.4+RX]; %you want to match these to the size of your imaging area
    slmYrange = [0.25 0.85];%9/21/20 [.05 0.9];%9/19/19 [.01 .7];% [0.075 0.85];%[0.5-RY 0.5+RY];
    
    % set Z range
    slmZrange = [-.09 .02];%3/11/21 IO.  %%[-0.02 0.07];%9/19/19 % [-0.1 0.15];
    % set to be approx -100 um (above SI) to + 22 (below SI) 3/10/21 WH
    % (new SLM)
    
    dummy = rand;
    
    slmCoords=zeros(4,npts);
    for i =1:npts
        slmCoords(:,i) = [...
            rand*(slmXrange(2)-slmXrange(1))+slmXrange(1),...
            rand*(slmYrange(2)-slmYrange(1))+slmYrange(1),...
            rand*(slmZrange(2)-slmZrange(1))+slmZrange(1),...
            1];
    end
    
    figure(1);scatter3(slmCoords(1,:),slmCoords(2,:),slmCoords(3,:),'o')
    drawnow;
    %%compile random holograms
    
    slmCoords(3,:) = round(slmCoords(3,:),3); %Added 3/15/21 by Ian for faster compute times

    
    
    disp('Compiling Holograms...')
    t = tic;
    try
        [ multiHolo,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoords' );
        multiPts = npts;
    catch
        multiPts = 100;%round(npts/2);
        disp('Could not create multi holo, trying with fewer points')
        [ multiHolo,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoords(:,1:multiPts));
    end
    fprintf(['Multi target Holo took ' num2str(toc(t)) 's\n'])
    multiCompileT = toc(t);
    
    % querry = input('do you want to check the range of the holos. turn off blasting then (1 yes, 0 no)');
    %
    % if querry ==1
    %      %%Check Range of multi holo
    %  disp('Starting with shutter closed, will display multiHolo. check that the range is appropriate on the basler');
    %   Function_Feed_SLM( Setup.SLM, multiHolo);
    %   mssend(masterSocket,[pwr/1000 1 multiPts]);
    %   function_BasPreview(Setup); %function_Basler_Preview(Setup, 5);
    % mssend(masterSocket,[0 1 1]);
    % end
    
   disp('Compiling Single Holograms')
%    disp('|------------------------------|')
%    fprintf('|')
%    tikmark = round(npts/30);
    parfor i =1:npts
        t=tic;
        fprintf(['Holo ' num2str(i)]);
        subcoordinates = slmCoords(:,i);
        
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        hololist(:,:,i)=Hologram;
        fprintf([' took ' num2str(toc(t)) 's\n']);
        
%         if mod(i,tikmark)==0
%         fprintf('.')
%         end
    end
%     fprintf('|\n')
    
    out.hololist = hololist;
    out.slmCoords = slmCoords;
    out.multiHolo = multiHolo;
    save('tempHololist4_will.mat','out');
else
    disp('Reloading old Holograms...')
    try
        load('tempHololist3.mat','out');
    catch
        [f, p] =uigetfile;
        load(fullfile(p,f),'out');
    end
    hololist = out.hololist;
    slmCoords = out.slmCoords;
    npts = size(slmCoords,2);
    multiHolo = out.multiHolo;
    figure(1);scatter3(slmCoords(1,:),slmCoords(2,:),slmCoords(3,:),'o')
    
end

disp(['Done compiling holograms. Took ' num2str(toc(tSingleCompile)) 's']);

singleCompileT = toc(tSingleCompile);

out.hololist=[];

%% Collect background frames for signal to noise testing
disp('Collecting Background Frames');

nBackgroundFrames = 20;

Bgdframe = function_BasGetFrame(Setup,nBackgroundFrames);% function_Basler_get_frames(Setup, nBackgroundFrames );
Bgd = uint8(mean(Bgdframe,3));
BGD = mean(Bgdframe,3);
meanBgd = mean(single(Bgdframe(:)));
stdBgd =  std(single(Bgdframe(:)));

threshHold = meanBgd+3*stdBgd;

fprintf(['3\x03c3 above mean threshold ' num2str(threshHold,4) '\n'])

%% Scan Image Planes Calibration
disp('Begining SI Depth calibration, we do this first incase spots burn holes with holograms')
tSI=tic;

zsToUse = linspace(0,70,11);% %70 was about 125um on 3/11/21 %Newer optotune has more normal ranges 9/28/29; New Optotune has different range 9/19/19; [0:10:89]; %Scan Image Maxes out at 89

SIUZ = -25:5:125;% linspace(-120,200,SIpts);
SIpts = numel(SIUZ);

%generate xy grid
%this is used bc SI depths are not necessarily parallel to camera depths
%but SI just generates a sheet of illumination
%therefore, we want to check a grid of spots that we will get depth info
%for
sz = size(Bgd);
gridpts = 25;
xs = round(linspace(1,sz(1),gridpts+2));
ys = round(linspace(1,sz(2),gridpts+2));

xs([1 end])=[];
ys([1 end])=[];
range =15;

%frames to average for image (orig 6) %added 7/15/2020 -Ian
framesToAcquire = 20;

clear dimx dimy XYSI
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

SIVals = zeros([SIpts c numel(zsToUse)]);

for k =1:numel(zsToUse)
    t=tic;
    z = zsToUse(k);
    fprintf(['Testing plane ' num2str(z) ': ']);
    
    
    mssend(SISocket,[z 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(SISocket,0.01);
    end
    
    
    dataUZ = zeros([sz SIpts]);
    for i = 1:numel(SIUZ)
        fprintf([num2str(round(SIUZ(i))) ' ']);
        
        currentPosition = getPosition(Sutter.obj);
        position = Sutter.Reference;
        position(3) = position(3) + (SIUZ(i));
        diff = currentPosition(3)-position(3);
        %           tic;
        moveTime=moveTo(Sutter.obj,position);
        %           toc;
        if i==1
            pause(1)
        else
            pause(0.1);
        end
        
        %change this part to change number of frames acquired
        % changed to 10 on 1/28/20 by WH to get better optotune calib
        % now a variable above... 7/15/2020 -Ian
        frame = function_BasGetFrame(Setup,framesToAcquire);%function_Basler_get_frames(Setup, 3 );
        %         frame = uint8(mean(frame,3));
        %         frame =  max(frame-Bgd,0);
        
        frame = (mean(frame,3)); %no cast to uint8
        frame =  max(frame-BGD,0);
        
        
        
        frame = imgaussfilt(frame,2); %removed because its binarizing for
%         some reason?
        dataUZ(:,:,i) =  frame;
        
        %          figure(1);
        %          subplot(1,2,1);
        %          imagesc(frame);
        %          subplot(1,2,2);
        %          imagesc(nanmean(dataUZ,3));
        %          drawnow
    end
    position = Sutter.Reference;
    moveTime=moveTo(Sutter.obj,position);
    pause(0.1)
    
    mssend(SISocket,[z 0]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(SISocket,0.01);
    end
    
    
    for i =1:c
%         temp = dataUZ(dimx(:,i),dimy(:,i),:);
%         SIVals(:,i,k) = mean(temp(:));
        SIVals(:,i,k) = squeeze(mean(mean(dataUZ(dimx(:,i),dimy(:,i),:))));
    end
    
    
    disp([' Took ' num2str(toc(t)) 's']);
    
    
end

mssend(SISocket,'end');

disp(['Scanimage calibration done whole thing took ' num2str(toc(tSI)) 's']);
siT=toc(tSI);

%%

disp('No Longer Putting off the actual analysis until later, Just saving for now')
out.SIVals =SIVals;
out.XYSI =XYSI;
out.zsToUse =zsToUse;
out.SIUZ = SIUZ;

save('TempSIAlign.mat','out')

%% First Fits
%Extract data to fit OptotuneZ as a function of camera XYZ

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
%%
SIpeakVal = b1;
SIpeakDepth = b2;

% thresholdModifier = 5; 1.5; %Ian add 9/13/19 orig 1.5
% excl = SIpeakVal<(threshHold/thresholdModifier);

SIThreshHold = 3*stdBgd/sqrt(nBackgroundFrames + framesToAcquire);
%  SIThreshHold = 0.3;
excl = SIpeakVal<SIThreshHold;%changed to reflect difference better.\
% excl = SIpeakVal<SIThreshHold;

disp([num2str(numel(SIpeakDepth)) ' points total before exclusions'])
disp([num2str(sum(excl(:))) ' points excluded b/c below threshold'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;

% excl = SIpeakVal>255;
% SIpeakVal(excl)=nan;
% SIpeakDepth(excl)=nan;

excl = SIpeakDepth<-50 | SIpeakDepth>150; %upper bound added 7/15/2020 -Ian

disp([num2str(sum(excl(:))) ' points excluded b/c too deep'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;
disp([num2str(sum(~isnan(SIpeakDepth), 'all')) ' points remaining'])


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
f2.Position = [0 0.5 0.5 0.45];
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
f3.Position = [0.5 0.5 0.5 0.45];
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

%% Coarse Data
% npts=200;

disp('Begining Coarse Holo spot finding')
coarsePts = 9; %odd number please
coarseUZ = linspace(-50,150,coarsePts);
mssend(masterSocket,[0 1 1]);

invar='flush';
while ~isempty(invar)
    invar = msrecv(masterSocket,0.01);
end

vals = nan(coarsePts,npts);
xyLoc = nan(2,npts);
tstart=tic;
sz = size(Bgd);
sizeFactor = 1;% changed bc camera size change 7/16/2020 by Ian %will have to manually test that this is scalable by 4
newSize = sz / sizeFactor;

% dataUZ2 = uint8(zeros([newSize  numel(coarseUZ) npts]));
dataUZ2 = zeros([newSize  numel(coarseUZ) npts], 'uint8');
maxProjections=uint8(zeros([newSize  npts]));

range=round(16 / sizeFactor);

numFramesCoarseHolo = 3; %number of frames to collect here. added 7/16/2020 -Ian

for i = 1:numel(coarseUZ)
    fprintf(['First Pass Holo, Depth: ' num2str(coarseUZ(i)) '. Holo : '])
    t = tic;
    
    currentPosition = getPosition(Sutter.obj);
    position = Sutter.Reference;
    position(3) = position(3) + (coarseUZ(i));
    diff = currentPosition(3)-position(3);
    moveTime=moveTo(Sutter.obj,position);
    
    if i==1
        pause(1)
    else
        pause(0.1);
    end
    
    for k=1:npts
        fprintf([num2str(k) ' ']);
        
        if mod(k,25)==0
            fprintf('\n')
        end
        
        
        Function_Feed_SLM( Setup.SLM, hololist(:,:,k));
        
        mssend(masterSocket,[pwr/1000 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        frame = function_BasGetFrame(Setup,numFramesCoarseHolo);% numFramesCoarseHolo added to be used elsewhere 7/16/2020 -Ian
        frame = uint8(mean(frame,3));
        
        mssend(masterSocket,[0 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        frame =  max(frame-Bgd,0);
        frame = imgaussfilt(frame,2);
        frame = imresize(frame,newSize);
        dataUZ2(:,:,i,k) =  frame;
        
        %          figure(1);
        %          subplot(1,2,1);
        %          imagesc(frame);
        %          subplot(1,2,2);
        %          imagesc(nanmean(dataUZ,3));
        %          drawnow
    end
    fprintf(['\nPlane Took ' num2str(toc(t)) ' seconds\n'])
    
end

position = Sutter.Reference;
moveTime=moveTo(Sutter.obj,position);
pause(0.1)
%%
disp('Calculating Depths and Vals')
range=ceil(15/sizeFactor);%Should be scaled by sizeFactor, also shrunk -Ian 7/16/2020 %Range for Hologram analysis window Changed to 5 9/16/19 by Ian
for k=1:npts
    dataUZ = dataUZ2(:,:,:,k);
    mxProj = max(dataUZ,[],3);
    [ x,y ] =function_findcenter(mxProj );
    xyLoc(:,k) = [x,y]*sizeFactor;
    
    maxProjections(:,:,k)=mxProj;
    
    dimx = max((x-range),1):min((x+range),size(mxProj,1));
    dimy =  max((y-range),1):min((y+range),size(mxProj,2));
    
    thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
    vals(:,k) = thisStack;
    depthIndex = find(thisStack == max(thisStack),1);
    
    fprintf(['Spot ' num2str(k) ' centered at depth ' num2str(round(coarseUZ(depthIndex)))...
        'um. Value: ' num2str(round(vals(depthIndex,k))) '\n']);
end
%     fprintf(['Took ' num2str(toc(t),2) 's\n']);

fprintf(['All Done. Total Took ' num2str(toc(tstart)) 's\n']);
coarseT = toc(tstart);

%% Second pass, multi-target version

% assign search params
finePts = 13; % odd number please
fineRange = 40;

disp('Begin multi-target z search...')
multi_time = tic;

% flush the socket
flushSocket(masterSocket)

% Generate multi-target holos based off coarse search data
targ_time = tic;


[coarseVal, coarseZidx] =max(vals,[],1);
zDepthVal = coarseUZ(coarseZidx);
zdepths = unique(zDepthVal);
n_planes = numel(zdepths);

coarseInclusionThreshold = 4*stdBgd/sqrt(numFramesCoarseHolo + nBackgroundFrames); %inclusion threshold added based on frames acquired; more stringent then SI. Added 7/16/2020 -Ian
zDepthVal(coarseVal<coarseInclusionThreshold)=NaN;

xyzLoc = [xyLoc;zDepthVal]; %fix this later (?)

%slmMultiCoords = nan(4, n_targs, n_holos);
%(4,:,:) = 1;  % 4th position is weight

clear slmMultiCoords basCoords targ_list targListIndiv slmMultiCoordsIndiv tempTargList

for i=1:n_planes % this will be the number of holograms
    % index the z depth
    z = zdepths(i);
    targ_idx = find(xyzLoc(3,:)==z);
    slmMultiCoords{i} = slmCoords(:,targ_idx);
    basCoords{i} = xyzLoc(:,targ_idx);
    targ_list{i} = targ_idx;
end

for i=1:n_planes
    % get real and slm coords from coarse
    dist = pdist2(basCoords{i}',basCoords{i}');
    %             dist(find(diag(diag(dist))))=NaN;
    temp =rand(size(dist,1));
    dist(find(diag(diag(temp))))=nan;  %#ok<FNDSB>
    tempTargList = 1:numel(targ_list{i});
    iterCounter =0;
    multiHoloCounter = 0;
    keepGoing=1;
    iterationsBeforeStop =1000;
    distanceThreshold = 30; %changed from 50 on 7/15/20 bc new cam 
    size_of_holo = 20;%changed from 25 on 7/15/20 bc new cam 
    doThisOnce =0;
    slmMultiCoordsIndiv{i} =[];
    targListIndiv{i}=[];
    
    while keepGoing
        iterCounter=iterCounter+1;
        if numel(tempTargList) <= size_of_holo
            testIdx = tempTargList;
            IdxofTempTargetList = 1:numel(tempTargList);
            keepGoing =0;
            %                 elseif numel(tempTargList)==2
            %                     return
        else
            IdxofTempTargetList = randperm(numel(tempTargList),size_of_holo);
            testIdx = tempTargList(IdxofTempTargetList);
        end
        
        %test if good
        subDist = dist(testIdx,testIdx);
        if any(subDist(:)<distanceThreshold)
            good =0;
        else
            good =1;
        end
        
        %complexish
        %                 a = find(any(subDist<distanceThreshold));
        %                 toBeKilled = a(2:end);
        %                 testIdx(toBeKilled) = [];
        %                 IdxofTempTargetList(toBeKilled) = [];
        
        if good
            multiHoloCounter=multiHoloCounter+1;
            slmMultiCoordsIndiv{i}{multiHoloCounter} = slmMultiCoords{i}(:,testIdx);
            targListIndiv{i}{multiHoloCounter} = targ_list{i}(testIdx) ;
            iterCounter=0;
            tempTargList(IdxofTempTargetList)=[];
        else
            if iterCounter>iterationsBeforeStop && doThisOnce
                keepGoing=0;
            elseif iterCounter>iterationsBeforeStop
                size_of_holo=max(round(size_of_holo/2),3);
                iterCounter=0;
                doThisOnce=1;
            end
        end 
    end
end


             

    % save to a struct for reference later
%     mh(i).slm = slmMultiCoords(:,:,i);
%     mh(i).real = basCoords;
%     mh(i).idx = holo_idx;
    
disp('Setting up stuff for multi-targets...');
% [Setup ] = function_loadparameters();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
Setup.useGPU =1;

cores=2;

if cores > 1
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=cores
        delete(p);
        parpool(cores);
    end
end

% make the holos
clear slmShootCoords
holo_time = tic;
disp('Compiling holograms...')
planes = numel(slmMultiCoordsIndiv);

for i=1:planes
    pt = tic;
    holos_this_plane = numel(slmMultiCoordsIndiv{i});
    
    parfor k=1:holos_this_plane
        ht = tic;
        [ mtholo, Reconstruction, Masksg ] = function_Make_3D_SHOT_Holos(Setup,slmMultiCoordsIndiv{i}{k}');
        
        mtholo_temp(k,:,:) = mtholo;
        %slm_temp(k,:,:) = slmMultiCoordsIndiv{i}{k}
        %disp(['Holo ' num2str(k) ' of ' num2str(holos_this_plane) ' done!  Took ' num2str(toc(ht)) 's'])
    end
    multiHolos{i} = mtholo_temp;
    disp(['Plane ' num2str(i) ' of ' num2str(planes) ' done!  Took ' num2str(toc(pt)) 's'])
end
disp(['Done. Took ' num2str(toc(holo_time)) 's'])

% multiHolos{plane}{holo, pixel, pixel}

disp(['took ' num2str(toc(targ_time)) 's to compile multi target holos']) 

out.hololist = hololist;
out.slmCoords = slmCoords;
out.multiHolo = multiHolo;
%save('tempHololist4_will_multi.mat','out');

%%
clear peakValue peakDepth peakFWHM
%%
background = function_BasGetFrame(Setup,20); % changed from 3 7/15/20
range = 6;
box_range = 20; % 7/15/20 changed from 50 to 20 distance threshold is set to 50, this must be less to avoid trying to fit 2 holos
disp('shootin!')
% for every  plane

for i = 1:planes%1:planes
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
        Function_Feed_SLM(Setup.SLM, multiHolos{i}(j,:,:));
        
        target_ref = targListIndiv{i}{j}(1);
        expected_z = xyzLoc(3,target_ref);
        %expected_xy = xyzLoc(1:2, target_ref);
        %expected_xyz = xyzLoc(:,target_ref);
        
        fineUZ = linspace(expected_z-fineRange,expected_z+fineRange,finePts);
        dataUZ = uint8(nan([size(Bgdframe(:,:,1))  finePts]));
        
        fprintf('Depth: ')
        
        % for every sutter z plane
        for k = 1:finePts
            fprintf([num2str(round(fineUZ(k))) ' ']);
            % move the sutter
            currentPosition = getPosition(Sutter.obj);
            position = Sutter.Reference;
            position(3) = position(3) + (fineUZ(k));
            diff = currentPosition(3)-position(3);
            moveTime=moveTo(Sutter.obj,position);

            if i==1
                pause(1)
            else
                pause(0.1);
            end
           
            requestPower(multi_pwr,masterSocket)
            
            % grab a frame, convert to uint8
            frame = function_BasGetFrame(Setup,10); %changed from 5 on 7/15/20
            frame = uint8(mean(frame,3));
            
            % turn off the laser
            requestPower(0,masterSocket)
            
            % subtract the background and filter
            frame =  max(frame-Bgd,0);
            frame = imgaussfilt(frame,2);
            % store into dataUZ(x,y,z-plane)
            dataUZ(:,:,k) =  frame;
            dataUZ3{i}{j} = dataUZ;
            fineUZ3{i}{j} = fineUZ;
        end
        
        % move sutter back to reference
        position = Sutter.Reference;
        moveTime=moveTo(Sutter.obj,position);
        pause(0.1)
        
        % OK, now parse the basler data in expected holo spots
        %targVals = nan(numel(fineZ),n_targs);
        for targ = 1:size(slmMultiCoordsIndiv{i}{j},2)
           
            target_ref = targListIndiv{i}{j}(targ);
            expected_xyz = xyzLoc(:,target_ref);
            [x, y] = size(Bgd);
            
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
                targ_stack = double(squeeze(max(max(dataUZ(targX,targY,:)))));
                mxProj = max(dataUZ(targX,targY,:),[],3);
                [ holo_x,holo_y ] =function_findcenter(mxProj );
                xyFine{i}{j}(:,targ) = [holo_x,holo_y];
            catch
                targ_stack = nan(finePts,1);
                xyFine{i}{j}(targ) = nan;
            end
           
%             try
%                 % method 2 - try to redected XY from multi-target hologram
%                 dimx = max((holo_x-range),1):min((holo_x+range),size(frame,1));
%                 dimy =  max((holo_y-range),1):min((holo_y+range),size(frame,2));
%                 targ_stack2 = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
%                 mxProj2 = max(dataUZ(dimX,dimY,:),[],3);
%                 [ holo_x2,holo_y2 ] =function_findcenter(mxProj2 );
%                 xyFine2{i}{j}(:,targ) = [holo_x2,holo_y2];
%                 
%             catch
%                 targ_stack2 = nan(finePts, 1);
%                 xyFine2{i}{j}(targ) = nan;
%             end
            
                % this is wrong
                %targVals(:,targ) = targ_stack;
                %this does not work yet
                %depthIndex = find(targ_stack == max(targVals),1);

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
            
%             try
%                 ff = fit(fineUZ', targ_stack2, 'gauss1');
%                 peakValue2{i}{j}(targ) = ff.a1;
%                 peakDepth2{i}{j}(targ) = ff.b1;
%                 peakFWHM2{i}{j}(targ) = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
%             catch
%                 disp(['Error on fit (method 2)! Holo: ', num2str(j), ' Target: ', num2str(targ)])
%                 peakValue2{i}{j}(targ) = NaN;
%                 peakDepth2{i}{j}(targ) = NaN;
%                 peakFWHM2{i}{j}(targ) = NaN;
%             end
        
        end
        
        fprintf('\n')
        disp(['Holo ' num2str(j) ' took ' num2str(toc(holo_time)) 's'])
    end
    
    disp(['Plane ' num2str(i) ' took ' num2str(toc(plane_time)) 's'])
end

fprintf(['All Done. Total Took ' num2str(toc(multi_time)) 's\n']);
multiT = toc(multi_time);

% %% reshape matrix of values into vectors
% peakValueList = peakValue(:);
% peakDepthList = peakDepth(:);
% peakFWHMList = peakFWHM(:);
% 
% [mx, mxi] = max(vals);
% % current implementation does not threshold, should be included but more
% % complicated here
% c=0;
% clear slmXYZ basXYZ basVal
% % for a = 1:npts %prob works but not always 750 pts bc of exlcusings
% %     c = c+1;
% for i=1:planes
%     holos_this_plane = numel(slmMultiCoordsIndiv{i});
%     for j=1:holos_this_plane
%         for targ = 1:size(slmMultiCoordsIndiv{i}{j},2)
%             c = c+1;
%             slmXYZ(:,c) = slmMultiCoordsIndiv{i}{j}(:,targ);
%             target_ref = targListIndiv{i}{j}(targ);
%             basXYZ2(1:2,c) = xyFine{i}{j}(targ);
%             basXYZ(1:2,c) = xyzLoc(1:2,target_ref);
%             basXYZ(3,c) = peakDepth{i}{j}(targ);
%             basVal(c) = peakValue{i}{j}(targ);
%         end
%     end
% end
%%
%% reshape matrix of values into vectors
peakValueList = peakValue(:);
peakDepthList = peakDepth(:);
peakFWHMList = peakFWHM(:);

[mx, mxi] = max(vals);
% current implementation does not threshold, should be included but more
% complicated here
c=0;
clear slmXYZ basXYZ1 basVal1
% for a = 1:npts %prob works but not always 750 pts bc of exlcusings
%     c = c+1;
for i=1:planes
    %     holos_this_plane = numel(slmMultiCoordsIndiv{i});
    
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
            
            % approach 1
%             basXYZ1_fine(1:2,c) = xyFine{i}{j}(targ);
            basXYZ1(1:2,c) = xyzLoc(1:2,target_ref);
            basXYZ1(3,c) = peakDepth{i}{j}(targ);
            basVal1(c) = peakValue{i}{j}(targ);
            
%             % approach 2
%             basXYZ2(1:2,c) = xyFine2{i}{j}(targ);
%             basXYZ2(3,c) = peakDepth2{i}{j}(targ);
%             basVal2(c) = peakValue2{i}{j}(targ);

FWHMval(c) = peakFWHM{i}{j}(targ); %added 7/10/20 =Ian
         end
    end
end

% figure(41)
% scatter3(basXYZ1(1,:),basXYZ1(2,:),basXYZ1(3,:), 'o', 'k')
% hold on
% scatter3(basXYZ2(1,:),basXYZ2(2,:),basXYZ2(3,:),'*','r')
% legend('Method 1', 'Method 2')
% title('Camera XYZ detected locations')
% 
% figure(45)
% scatter3(basXYZ1_fine(1,:),basXYZ1_fine(2,:),basXYZ1(3,:), 'filled', 'r', 'MarkerFaceAlpha',0.7)
% hold on
% scatter3(basXYZ1(1,:),basXYZ1(2,:),basXYZ1(3,:), 'filled', 'g', 'MarkerFaceAlpha',0.7)
% legend('XY from fine', 'XY from coarse')
% title('Camera XYZ detected locations, same method')


%% Choose your favorite method of getting XYZ coords
approach = 1;
disp(['you chose approach ' num2str(approach)])
switch approach
    case 1
        basXYZ = basXYZ1;
        basVal = basVal1;
    case 2
        basXYZ = basXYZ2;
        basVal = basVal2;
end


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
FWHMVal = FWHMBackup2;%added 7/20/2020 -Ian

excludeTrials = all(basXYZ(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

excludeTrials = excludeTrials | basVal>260; %max of this camera is 255

basDimensions = size(Bgdframe);
excludeTrials = excludeTrials | basXYZ(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ(3,:)<-50; %9/19/19 Ian Added to remove systematic low fits
excludeTrials = excludeTrials | basXYZ(3,:)>190;


excludeTrials = excludeTrials | any(isnan(basXYZ(:,:)));
excludeTrials = excludeTrials | basVal<1; %8/3 hayley add 5; Ian ammend to 1 9/13
excludeTrials = excludeTrials | basVal>(mean(basVal)+3*std(basVal)); %9/13/19 Ian Add

slmXYZBackup = slmXYZ(:,~excludeTrials);
basXYZBackup = basXYZ(:,~excludeTrials);
basValBackup = basVal(:,~excludeTrials);
FWHMValBackup = FWHMVal(~excludeTrials); % added 7/20/2020 -Ian
%basValBackup = basValBackup(:,1:386); % WH add to exlude trials with water loss
%slmXYZBackup = slmXYZBackup(:,1:386); % did this on 1/29/20
%basXYZBackup = basXYZBackup(:,1:386);
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
plot(basValBackup,'o')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')






%% fit SLM to Camera
%use model terms

basXYZ = basXYZBackup;
slmXYZ = slmXYZBackup;
basVal = basValBackup;

disp('Fitting SLM to Camera')
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1 ; 1 1 1 ;...
    2 0 0; 0 2 0; 0 0 2;  ...
    2 0 1; 2 1 0; 0 2 1; 1 2 0; 0 1 2;  1 0 2; ... ];  %XY spatial calibration model for Power interpolations
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2;];
reOrder = randperm(size(slmXYZ,2));
slmXYZ = slmXYZ(:,reOrder);
basXYZ = basXYZ(:,reOrder);

holdback = 100;%50;

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = (basXYZ(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

errScalar = 2.5; %2.8;%2.5;
figure(1286);clf;subplot(1,2,1)
[SLMtoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
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

% pxPerMu = size(frame,1) / 1000; %really rough approximate of imaging size
pxPerMu = 0.57723;

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

%  camToSLM = function_3DCoC(refAsk,refGet,modelterms);

subplot(1,2,2)
[camToSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
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

% RMS = sqrt(sum((Get-True).^2,2));
% meanRMS = nanmean(RMS);
% 
% disp(['The RMS error: ' num2str(meanRMS) ' SLM units for Camera to SLM']);




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
%  modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%      1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%      2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;   ];  %XY spatial calibration model for Power interpolations
slmXYZ = slmXYZBackup;
basVal = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

intVal = basVal;
intVal = sqrt(intVal); %convert fluorescence intensity (2P) to 1P illumination intensity
intVal=intVal./max(intVal(:));

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = intVal(1:end-holdback);

SLMtoPower =  polyfitn(refAsk,refGet,modelterms);

Ask = (slmXYZ(1:3,end-holdback:end))';
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

%% Now using these CoC lets create holograms that shoot a pattern into a field of view
disp('Picking Holes to Burn')


%module to restrict burn grid to most likely SI FOV.
%added 7/20/19 -Ian
binarySI = ~isnan(SIpeakVal);
SImatchProb = mean(binarySI');%probability that point was detected aka that was in SI range

SImatchXY = camXYZ(1:2,1:625); %location of points in XYZ

figure(8);clf;
s=scatter(SImatchXY(1,:),SImatchXY(2,:),[],SImatchProb,'filled');

SImatchThreshold = 0.5; % threshold for being in SI FOV (set to 0 to take whole range)

SIx = SImatchXY(1,SImatchProb>SImatchThreshold);
SIy = SImatchXY(2,SImatchProb>SImatchThreshold);

SIboundary = boundary(SIx',SIy');
hold on
p=plot(SIx(SIboundary),SIy(SIboundary));

sz = size(Bgd);
%for fine targets
bufferMargin = 0.05; %fraction of total area as buffer
SImatchRangeXforFine = [max(min(SIx(SIboundary))-sz(1)*bufferMargin,1) ...
    min(max(SIx(SIboundary))+sz(1)*bufferMargin,sz(1))];

SImatchRangeYforFine = [max(min(SIy(SIboundary))-sz(2)*bufferMargin,1) ...
    min(max(SIy(SIboundary))+sz(2)*bufferMargin,sz(2))];
r = rectangle('position',...
    [SImatchRangeXforFine(1) SImatchRangeYforFine(1) SImatchRangeXforFine(2)-SImatchRangeXforFine(1) SImatchRangeYforFine(2)-SImatchRangeYforFine(1)]);
r.EdgeColor='g';

%for hole burning
bufferMargin = 0.02; %fraction of total area as buffer
SImatchRangeX = [max(min(SIx(SIboundary))-sz(1)*bufferMargin,1) ...
    min(max(SIx(SIboundary))+sz(1)*bufferMargin,sz(1))];

SImatchRangeY = [max(min(SIy(SIboundary))-sz(2)*bufferMargin,1) ...
    min(max(SIy(SIboundary))+sz(2)*bufferMargin,sz(2))];

r = rectangle('position',...
    [SImatchRangeX(1) SImatchRangeY(1) SImatchRangeX(2)-SImatchRangeX(1) SImatchRangeY(2)-SImatchRangeY(1)]);
r.EdgeColor='r';

rline = line(NaN,NaN,'LineWidth',1','LineStyle', '-','color','r');
gline = line(NaN,NaN,'LineWidth',1','LineStyle', '-','color','g');

colorbar
legend('Prob of SI FOV','Detected SI FOV','Burn Boundary Box','Calib Boundary Box')
xlabel('Camera X pixels')
ylabel('Camera Y Pixels')

nBurnGrid = 8; %number of points in the burn grid
xpts = linspace(SImatchRangeX(1),SImatchRangeX(2),nBurnGrid);
ypts = linspace(SImatchRangeY(1),SImatchRangeY(2),nBurnGrid);

XYpts =[];
for i=1:nBurnGrid
    for k=1:nBurnGrid
        XYpts(:,end+1) = [xpts(i) ypts(k)];
    end
end

zsToBlast = zsToUse;% match to SI Calib %linspace(0,90,11);% Changed to account for newer optotune 9/28/20; Changed to account for new optotune range 9/19/19 by Ian 0:10:80; %OptoPlanes to Blast
interXdist = xpts(2)-xpts(1);
%  xOff = round(interXdist/numel(zsToBlast));
interYdist = ypts(2)-ypts(1);
%  yOff = round(interYdist/numel(zsToBlast));

gridSide = ceil(sqrt(numel(zsToBlast)));
xOff = round(interXdist/gridSide);
yOff = round(interYdist/gridSide);

%Turn into a more unique looking pattern
numPts = size(XYpts,2);
FractionOmit = 0.1; %changed down to 10% from 25% bc not really needed. 9/19/19 by Ian
XYpts(:,randperm(numPts,round(numPts*FractionOmit)))=[];
XYpts = reshape(XYpts,[2 numel(XYpts)/2]);

disp([num2str(size(XYpts,2)) ' points per plane selected. ' num2str(size(XYpts,2)*numel(zsToBlast)) ' total'])

%% Simulate and create new Fine POints
% do a CoC to get more points to shoot

denseFineTimer = tic; 

nSimulatedTargs = 10000;

% get basler targets to shoot in from basler range
% for X
a = min(basXYZBackup(1,:));
b = max(basXYZBackup(1,:));
% a =SImatchRangeXforFine(1);
% b = SImatchRangeXforFine(2);
r = (b-a).*rand(nSimulatedTargs,1) + a;
rX = round(r);

% for Y
a = min(basXYZBackup(2,:));
b = max(basXYZBackup(2,:));
r = (b-a).*rand(nSimulatedTargs,1) + a;
rY = round(r);

% for Z
a = min(basXYZBackup(3,:));
b = max(basXYZBackup(3,:));
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
multiholosize=20;
planes = 7; 
holosperplane = 5;
ntotalPoints = multiholosize * planes * holosperplane;
disp(['Using ' num2str(ntotalPoints) ' points in round 2.'])

slm_coords = {};
bas_coords = {};
%c = 0;

[idx, cent] = kmeans(expectBas(:,3), planes);
% figure
% scatter3(expBas(:,1), expBas(:,2), expBas(:,3), [], categorical(idx))

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

c = 0;
for i=1:length(slm_coords)
    for j=1:length(slm_coords{i})
        c = c + 1;
        ht = tic;
        disp(['Compiling multi-target hologram ' num2str(c)])
        
        thisCoord = slm_coords{i}{j};
        thisCoord(:,3) = round(thisCoord(:,3),3); %Added 3/15/21 by Ian for faster compute times 

        [ mtholo, Reconstruction, Masksg ] = function_Make_3D_SHOT_Holos(Setup,thisCoord);
        holos2shoot{i}{j} = mtholo;
        disp(['done in ' num2str(toc(ht)) 's.'])
    end          
end

disp('now to shooting...') 

%% now repeat multi-target search with new holograms
clear peakValue4 peakDepth4 peakFWHM4 dataUZ4

%%

background = function_BasGetFrame(Setup,20);
range = 6;
box_range = 20; % distance threshold is set to 100
disp('shooting!')

planes = numel(slm_coords);
for i = 1:planes
    holos_this_plane = numel(slm_coords{i});
    disp(['Plane ' num2str(i) ' of ' num2str(planes)])
    
    % find the mean z of the holo targets and set a range around it
    % rearanging so it doesn't move unnescessarily
    meanz = mean(cellfun(@(x) mean(x(:,3)),bas_coords{i}));
%         minz = min(cellfun(@(x) min(x(:,3)),bas_coords{i}));
%         maxz = min(cellfun(@(x) max(x(:,3)),bas_coords{i}));

%      meanz = mean(bas_coords{i}{j}(:,3));
    fineUZ = linspace(meanz-fineRange, meanz+fineRange, finePts);

    % for every holo on that plane
%     for j = 1:holos_this_plane
        dataUZPlane = uint8(nan([size(Bgdframe(:,:,1)) finePts holos_this_plane]));
        
        % set power
%         multi_pwr = size(slm_coords{i}{j},1) * pwr;        
%         Function_Feed_SLM(Setup.SLM, holos2shoot{i}{j});
%         
        
        
        % for every sutter z plane
        fprintf('Depth: ')
        figure(4);clf;
        for k = 1:finePts
            fprintf([num2str(round(fineUZ(k))) ' ']);
            
            % move the sutter
            currentPosition = getPosition(Sutter.obj);
            position = Sutter.Reference;
            position(3) = position(3) + (fineUZ(k));
            diff = currentPosition(3)-position(3);
            moveTime=moveTo(Sutter.obj,position);

            if i==1
                pause(1)
            else
                pause(0.1);
            end
           
            a = floor(sqrt(holos_this_plane));
            b = ceil(holos_this_plane/a);
            ro = min([a b]);
            co = max([a b]); 
            

            for j= 1:holos_this_plane
                multi_pwr = size(slm_coords{i}{j},1) * pwr;
                Function_Feed_SLM(Setup.SLM, holos2shoot{i}{j});
                
                requestPower(multi_pwr,masterSocket)
                
                % grab a frame, convert to uint8
                frame = function_BasGetFrame(Setup,3);
                frame = uint8(mean(frame,3));
                
                % turn off the laser
                requestPower(0,masterSocket)
                
                % subtract the background and filter
                frame =  max(frame-Bgd,0);
                frame = imgaussfilt(frame,2);
                % store into dataUZ(x,y,z-plane)
                dataUZPlane(:,:,k,j) =  frame;
                %             dataUZ4{i}{j} = dataUZ;
                
                figure(4);
                subplot(ro,co,j)
                imagesc(frame)
                colorbar
                caxis([0 150]);
                title({['Live Data. Depth ' num2str(round(fineUZ(k)))] ; ['Plane: ' num2str(i) '. Set ' num2str(j)]})
                drawnow
       
            end
        end
        
        for j= 1:holos_this_plane           
        fineUZ4{i}{j} = fineUZ;
        dataUZ4{i}{j} = dataUZPlane(:,:,:,j); %brought out of for loop
        
        dataUZ = dataUZPlane(:,:,:,j);
        
        % move sutter back to reference
        position = Sutter.Reference;
        moveTime=moveTo(Sutter.obj,position);
        pause(0.1)
        
        % target parsing, might do later instead
        for targ = 1:size(slm_coords{i}{j},1)
            
            expected_xyz = bas_coords{i}{j}(targ,:);
            [x, y] = size(Bgd);
            
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
                targ_stack = double(squeeze(max(max(dataUZ(targX,targY,:)))));
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

excludeTrials = excludeTrials | basVal4>260; %max of this camera is 255

basDimensions = size(Bgdframe);
excludeTrials = excludeTrials | basXYZ4(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ4(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ4(3,:)<-50; %9/19/19 Ian Added to remove systematic low fits
excludeTrials = excludeTrials | basXYZ4(3,:)>190;


excludeTrials = excludeTrials | any(isnan(basXYZ4(:,:)));
excludeTrials = excludeTrials | basVal4<1; %8/3 hayley add 5; Ian ammend to 1 9/13
excludeTrials = excludeTrials | basVal4>(mean(basVal4)+3*std(basVal4)); %9/13/19 Ian Add

slmXYZBackup = slmXYZ4(:,~excludeTrials);
basXYZBackup = basXYZ4(:,~excludeTrials);
basValBackup = basVal4(:,~excludeTrials);
FWHMValBackup = FWHMVal4(~excludeTrials); % added 7/20/2020 -Ian
%basValBackup = basValBackup(:,1:386); % WH add to exlude trials with water loss
%slmXYZBackup = slmXYZBackup(:,1:386); % did this on 1/29/20
%basXYZBackup = basXYZBackup(:,1:386);
%%
f41=figure(1922);
f41.Units = 'Normalized';
f41.Position = [0.05 0.4 0.5 0.5];
subplot(1,2,1)
% scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 50, 'k', 'Filled', 'MarkerFaceAlpha',0.5)
% hold on
% scatter3(xyzLoc(1,:), xyzLoc(2,:), xyzLoc(3,:), 70,  'r', 'Filled', 'MarkerFaceAlpha', 0.7)
% legend('Fine','Coarse')
% title('Detected basXYZs')
% subplot(1,2,2)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 75, basValBackup, 'Filled')
colorbar
colormap default
title({'Second Denser Fine';'basXYZ and basVals (fine)'})

% f42=figure(42);
% clf(42)
subplot(1,2,2)
plot(basValBackup,'o')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')






%% fit SLM to Camera
%use model terms

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

holdback = 100;%50;

refAsk = (slmXYZ4(1:3,1:end-holdback))';
refGet = (basXYZ4(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

errScalar = 2.5; %2.8;%2.5;
figure(1977);clf;subplot(1,2,1)
[SLMtoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
title('SLM to Cam v2')

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


Ask = (slmXYZ4(1:3,end-holdback:end))';
True = (basXYZ4(1:3,end-holdback:end))';
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

% pxPerMu = size(frame,1) / 1000; %really rough approximate of imaging size
pxPerMu = 0.57723;

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

subplot(1,2,2)
[camToSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
title('Cam to SLM v2')

Ask = (basXYZ4(1:3,end-holdback:end))';
True = (slmXYZ4(1:3,end-holdback:end))';
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

% RMS = sqrt(sum((Get-True).^2,2));
% meanRMS = nanmean(RMS);
% 
% disp(['The RMS error: ' num2str(meanRMS) ' SLM units for Camera to SLM']);




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
%  modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%      1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%      2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;   ];  %XY spatial calibration model for Power interpolations
slmXYZ4 = slmXYZBackup;
basVal4 = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

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

figure(1);clf
subplot(2,3,1)
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],intVal,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Measured Power (converted to 1p)')
colorbar
axis square

subplot(2,3,2)
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],polyvaln(SLMtoPower,slmXYZ4(1:3,:)'),'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated Power Norm.')
colorbar
axis square

subplot(2,3,4)
% scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal','filled');
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],basVal4,'filled');

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Raw Fluorescence')
colorbar
axis square

subplot(2,3,5)
c = sqrt((polyvaln(SLMtoPower,slmXYZ4(1:3,:)')-intVal').^2);
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error RMS (A.U.)')
colorbar
axis square

subplot(2,3,3)
c = (polyvaln(SLMtoPower,slmXYZ4(1:3,:)').^2);
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated 2P Power')
colorbar
axis square

subplot(2,3,6)
normVal = basVal4./max(basVal4(:));

c = (polyvaln(SLMtoPower,slmXYZ4(1:3,:)').^2)-normVal';
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
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
ylim([-75 150])
xlim([7.5 35])

refline(0,0)
refline(0,100)


subplot(1,2,2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],FWHM,'filled')
caxis([10 30])
h= colorbar;
xlabel('SLM X')
ylabel('SLM Y')
zlabel('SLM Z')
set(get(h,'label'),'string','FWHM \mum')

fprintf(['FWHM in the typical useable volume (0 to 100um) is: ' num2str(mean(FWHM(depth>0 & depth<100))) 'um\n'])



%% Plot Hole Burn Stuff
figure(6);
scatter(XYpts(1,:),XYpts(2,:),'o');


figure(4); clf


clear XYtarg SLMtarg
for i = 1:numel(zsToBlast)
    
    %offset the xy points each Z plane
    
    %      XYuse = bsxfun(@plus,XYpts,([xOff yOff].*(i-1))');
    
    a = mod(i-1,gridSide);
    b = floor((i-1)/gridSide);
    %      [a b]
    
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
    ExcludeBurns = estPower<0; %don't shoot if you don't have the power
    estSLM(ExcludeBurns,:)=[];
    estPower(ExcludeBurns)=[];
    XYuse(:,ExcludeBurns)=[];
    zOptoPlane(ExcludeBurns)=[];
    estCamZ(ExcludeBurns)=[];
%     
    
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
tCompileBurn = tic;

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

%%
disp('Blasting Holes for SI to SLM alignment, this will take about an hour and take 25Gb of space')
tBurn = tic;

%confirm that SI computer in eval Mode
mssend(SISocket,'1+2');
invar=[];
while ~strcmp(num2str(invar),'3') %stupid cludge so that [] read as false
    invar = msrecv(SISocket,0.01);
end
disp('linked')

%mssend(SISocket,'hSI.startGrab()');
%setup acquisition

numVol =5; %number of SI volumes to average
baseName = '''calib''';

mssend(SISocket,['hSI.hStackManager.arbitraryZs = [' num2str(zsToBlast) '];']);
mssend(SISocket,['hSI.hStackManager.numVolumes = [' num2str(numVol) '];']);
mssend(SISocket,'hSI.hStackManager.enable = 1 ;');

mssend(SISocket,'hSI.hBeams.pzAdjust = 0;');
mssend(SISocket,'hSI.hBeams.powers = 16;'); %power on SI laser. important no to use too much don't want to bleach

mssend(SISocket,'hSI.extTrigEnable = 0;'); %savign
mssend(SISocket,'hSI.hChannels.loggingEnable = 1;'); %savign
mssend(SISocket,'hSI.hScan2D.logFilePath = ''D:\Calib\Temp'';');
% mssend(SISocket,'hSI.hScan2D.logFileCounter = 1;');
mssend(SISocket,['hSI.hScan2D.logFileStem = ' baseName ';']);
mssend(SISocket,'hSI.hScan2D.logFileCounter = 1;');


% mssend(SISocket,'1+2');

mssend(SISocket,['hSICtl.updateView;']);

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
%clear invar
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
        disp(['Ready for Next'])
    else
        %             disp(invar)
    end
end


burnPowerMultiplier = 10;%change to 10 3/11/21 %previously 5; added by Ian 9/20/19
burnTime = 0.5; %in seconds, very rough and not precise

disp('Now Burning')

for k=1:numel(zsToBlast)%1:numel(zsToBlast)
    
    offset = round(meanCamZ(k));
    currentPosition = getPosition(Sutter.obj);
    position = Sutter.Reference;
    position(3) = position(3) + (offset);
    diff = currentPosition(3)-position(3);
    %           tic;
    moveTime=moveTo(Sutter.obj,position);
    %           toc;
    if k==1
        pause(1)
    else
        pause(0.1);
    end
    
    tempHololist=holos{k};
    
    for i=1:size(XYtarg{k},2)%1:size(XYuse,2)
        t=tic;
        fprintf(['Blasting Hole ' num2str(i) '. Depth ' num2str(zsToBlast(k))]);
        Function_Feed_SLM( Setup.SLM, tempHololist(:,:,i));
        
        DE = Diffraction{k}(i);
        if DE<.05 %if Diffraction efficiency too low just don't even burn %Ian 9/20/19
            DE=inf;
        end
        blastPower = pwr*burnPowerMultiplier /1000 /DE;
        
        if blastPower>2 %cap for errors, now using a high divided mode so might be high 
            blastPower =2;
        end
        
        stimT=tic;
        mssend(masterSocket,[blastPower 1 1]);
        while toc(stimT)<burnTime
        end
        mssend(masterSocket,[0 1 1]);
        
        %flush masterSocket %flush and handshake added 9/20/19 by Ian
        invar='flush';
        while ~isempty(invar)
            invar = msrecv(masterSocket,0.01);
        end
        %re send 0
        mssend(masterSocket,[0 1 1]);
        %check for handshake
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end

        
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

position = Sutter.Reference;
moveTime=moveTo(Sutter.obj,position);

burnT = toc(tBurn);
disp(['Done Burning. Took ' num2str(burnT) 's']);

disp('Done with Lasers and ScanImage now, you can turn it off')

%% Move file to modulation

disp('Moving files')
tMov = tic;

%on ScanImage Computer
% destination = '''F:\frankenshare\FrankenscopeCalib''' ;
destination = '''F:\Calib''';
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


%%
mssend(SISocket,'end');
%%

tLoad = tic;
% pth = 'W:\frankenshare\FrankenscopeCalib'; %On this compute
% pth ='F:\Temp'; %temp fix b/c frankenshare down
pth = 'F:\Calib';
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
            maskFR = imgaussfilt(Frame,3) - imgaussfilt(Frame,16);
            mask = maskFR > mean(maskFR(:))+6*std(maskFR(:));
            
            %remove the low frequency slide illumination differences
            filtNum = 4;
            frameFilt = imgaussfilt(Frame,filtNum);
            baseFilt = imgaussfilt(baseFrame,filtNum);
            
            
            toCalc = (baseFrame-baseFilt) - (Frame-frameFilt);
            toCalc(mask)=0;
            
%             testFr = Frames{k}(:,:,c-1) - Frame;
            [ x,y ] =function_findcenter(toCalc);
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
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
    
%     
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations


cam3XYZ = [XYtarg{:};];
SIXYZ = SIXYZbackup;

cam3XYZ=cam3XYZ(:,1:size(SIXYZ,2));

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;
cam3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];

% testSet = randperm(size(SIXYZ,2),25);
% otherSet = ones([size(SIXYZ,2) 1]);
% otherSet(testSet)=0;
% otherSet = logical(otherSet);

refAsk = SIXYZ(1:3,:)';
refGet = (cam3XYZ(1:3,:))';
errScalar = 2.;2.5;2.6;

figure(2594);clf
subplot(1,2,1)
[SItoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
title('SI to Cam')

subplot(1,2,2)
[CamToSI, trialN] = function_3DCoCIterative(refGet,refAsk,modelterms,errScalar,0);
title('Cam to SI')

CoC.CamToSI = CamToSI;
CoC.SItoCam = SItoCam;
out.CoC=CoC;

%% alternate calculation
% modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%     1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

tempSLM = cellfun(@(x) x',SLMtarg,'UniformOutput',false);
slm3XYZ = [tempSLM{:}];
SIXYZ = SIXYZbackup;

slm3XYZ=slm3XYZ(1:3,1:size(SIXYZ,2));

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;

slm3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];

refAsk = SIXYZ(1:3,:)';
refGet = (slm3XYZ(1:3,:))';
errScalar =2;2.5;2.5;

figure(2616);clf
subplot(1,2,1)
[SItoSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
title('SI to SLM')
subplot(1,2,2)
[SLMtoSI, trialN] = function_3DCoCIterative(refGet,refAsk,modelterms,errScalar,0);
title('SLM to SI')

CoC.SItoSLM = SItoSLM;
CoC.SLMtoSI = SLMtoSI;


%% Calculate round trip errors
numTest = 10000;

rangeX = [0 511];%[0 511];
rangeY = [0 511];%[0 511];
rangeZ = [0 70];% Make Sure to match this to the correct range for this optotune;

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

%% Save Output Function
pathToUse = 'C:\Users\Holography\Desktop\SLM_Management\Calib_Data';
disp('Saving...')
tSave = tic;

save(fullfile(pathToUse,[date '_Calib.mat']),'CoC')
save(fullfile(pathToUse,'ActiveCalib.mat'),'CoC')



times.saveT = toc(tSave);
times.loadT = loadT;
times.MovT = MovT;
times.burnT = burnT;
times.compileBurnT = compileBurnT;
times.fitsT = fitsT;
times.fineT = fineT; %Second Dense Fine
times.multiT = multiT; %first pass fine fit
times.coarseT = coarseT; %Coarse Fit
times.siT = siT;
times.singleCompileT = singleCompileT;
times.multiCompileT = multiCompileT;

totT = toc(tBegin);
times.totT = totT;
save(fullfile(pathToUse,'CalibWorkspace.mat'));
% save(fullfile(pathToUse,['CalibWorkspace_will_' date '.mat']))
disp(['Saving took ' num2str(toc(tSave)) 's']);

disp(['All Done, total time from begining was ' num2str(toc(tBegin)) 's. Bye!']);

%% use this to run
%[SLMXYZP] = function_SItoSLM(SIXYZ,CoC);




[Setup.SLM ] = Function_Stop_SLM( Setup.SLM );

try; function_close_sutter( Sutter ); end
try function_stopBasCam(Setup); end



        
