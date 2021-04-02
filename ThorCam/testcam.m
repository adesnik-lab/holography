
try stopCamera(tlCamera,other); end

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

addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\ThorCam\'));

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

try stopCamera(tlCamera,other); end
% [cam, other] = thorSetupCam;
cam = thorSetupCam;

Setup.cam = cam;

disp('Ready')
%% look for objective in 1p or you know... have it already set up
thorPreview(Setup);

%% Make mSocketConnections with DAQ and SI Computers

disp('Waiting for msocket communication From DAQ')
%then wait for a handshake
srvsock = mslisten(4205);
masterSocket = msaccept(srvsock,30);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);
%MasterIP = '128.32.177.217';
%masterSocket = msconnect(MasterIP,3002);

invar = [];

while ~strcmp(invar,'B');
    invar = msrecv(masterSocket,.5);
end;
disp('communication from Master To Holo Established');

%% Set Power Levels

pwr =40;25;%;14.3; 12.5;12; %40; %13 at full; 50 at 15 divided
disp(['individual hologram power set to ' num2str(pwr) 'mW']);

%%
disp('Find the spot and check if this is the right amount of power')
slmCoordsTemp = [0.5 0.4 0 1];%[0.276 .4633 .01676 1];%[0.4 0.75 0.01 1];
DEestimateTemp = DEfromSLMCoords(slmCoordsTemp); %
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimateTemp)])

[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoordsTemp );

blankHolo = zeros([1920 1152]);
%Function_Feed_SLM( Setup.SLM, blankHolo);

Function_Feed_SLM( Setup.SLM, HoloTemp);

mssend(masterSocket,[pwr/1000 1 1]);

thorPreview(cam);
% function_BasPreview(Setup);
mssend(masterSocket,[0 1 1]);

%% Collect background frames for signal to noise testing
disp('Collecting Background Frames');

nBackgroundFrames = 5;

Bgdframe = thorAcquireSimple(nBackgroundFrames, Setup);% function_Basler_get_frames(Setup, nBackgroundFrames );
Bgd = uint16(mean(Bgdframe,3));
meanBgd = mean(single(Bgdframe(:)));
stdBgd =  std(single(Bgdframe(:)));

threshHold = meanBgd+3*stdBgd;

fprintf(['3\x03c3 above mean threshold ' num2str(threshHold,4) '\n'])
%%
requestPower(40,masterSocket);

[ave_im] = thorAcquireSimple(5,tlCamera);
%%
clear aves mxs t aves2
% exposure_times = [10000, 100000, 200000, 500000];
exposure_times = [1e3, 5e3, 10e3, 50e3, 100e3, 500e3];
frames_avg = 1;

nBackgroundFrames = 5;
Bgdframe = thorAcquireSimple(nBackgroundFrames, 100e3, cam);% function_Basler_get_frames(Setup, nBackgroundFrames );
Bgd = mean(Bgdframe,3);

% 
% for i=1:numel(exposure_times)
%     requestPower(0,masterSocket)
%     requestPower(40,masterSocket);
%     tic
%     [ave_im, mx] = thorAcquireSimple(10,exposure_times(i), tlCamera);
%     t(i) = toc;
%     requestPower(0,masterSocket)
%     aves(i,:,:) = ave_im;
%     mxs(i) = mx;
% end

for i=1:numel(exposure_times)
    requestPower(0,masterSocket)
    requestPower(40,masterSocket);
    tic
    [frame, mx] = thorAcquireSimple(frames_avg ,exposure_times(i), cam);
    t(i) = toc;
    requestPower(0,masterSocket)
    frame =  max((frame)-Bgd,0);
    frame = imgaussfilt(frame,2);
    aves(i,:,:) = frame;
%     mxs(i) = mx;
    [x,y] = function_findcenter(frame);
    
    try
        aves2(i,:,:) = frame(x-20:x+20, y-20:y+20);
    catch
        aves2(i,:,:) = nan(41,41);
    end
    
    mxs(i) = max(aves2(i,:,:), [], 'all');
    
end

bits = (2^16)-1;

figure(920)
clf
yyaxis left
plot(exposure_times/1000, mxs,  '-o')
% plot(frames_avg, mxs/bits*100, '-o')
ylabel('px val')
hold on
yyaxis right
plot(exposure_times/1000, t*1000,  '-o')
% plot(frames_avg, t*1000,  '-o')
ylabel('time (ms)')
xlabel('exposure time (ms)')
set(gca, 'XScale', 'log')
hold off
name = ['plot_' num2str(frames_avg) '_frames.png'];
saveas(gcf, ['./thorcamtests/' name])


figure(930)
clf
nplts = numel(mxs);
c=0;
for i = 1:numel(mxs)
    c=c+1;
    subplot(2,nplts,c)
    imagesc(squeeze(aves(i,:,:)))
    colorbar
end
hold on
for i = 1:numel(mxs)
    c=c+1;
    subplot(2,nplts,c)
    imagesc(squeeze(aves2(i,:,:)))
    colorbar
end
hold off
name = ['holos_' num2str(frames_avg) '_frames.png'];
saveas(gcf, ['./thorcamtests/' name])

%%
%% Collect background frames for signal to noise testing
disp('Collecting Background Frames');

nBackgroundFrames = 5;

Bgdframe = thorAcquireSimple(nBackgroundFrames, tlCamera);% function_Basler_get_frames(Setup, nBackgroundFrames );
Bgd = uint16(mean(Bgdframe,3));
meanBgd = mean(single(Bgdframe(:)));
stdBgd =  std(single(Bgdframe(:)));

threshHold = meanBgd+3*stdBgd;

fprintf(['3\x03c3 above mean threshold ' num2str(threshHold,4) '\n'])

requestPower(40,masterSocket);

[ave_im] = thorAcquireSimple(5,tlCamera);
%%
clear aves mxs t aves2
% exposure_times = [10000, 100000, 200000, 500000];
exposure_times = [1e3, 5e3, 10e3, 50e3, 100e3, 500e3];
frames_avg = [1, 3, 5, 10];

nBackgroundFrames = 5;
Bgdframe = function_BasGetFrame(Setup, nBackgroundFrames);% function_Basler_get_frames(Setup, nBackgroundFrames );
Bgd = uint8(mean(Bgdframe,3));

% 
% for i=1:numel(exposure_times)
%     requestPower(0,masterSocket)
%     requestPower(40,masterSocket);
%     tic
%     [ave_im, mx] = thorAcquireSimple(10,exposure_times(i), tlCamera);
%     t(i) = toc;
%     requestPower(0,masterSocket)
%     aves(i,:,:) = ave_im;
%     mxs(i) = mx;
% end

for i=1:numel(frames_avg)
    requestPower(0,masterSocket)
    requestPower(40,masterSocket);
    [frame, tm] = function_BasGetFrame( Setup, frames_avg(i));
    t(i) = tm;
    requestPower(0,masterSocket)
    frame =  max(frame(:,:,1)-Bgd,0);
    frame = imgaussfilt(frame,2);
    aves(i,:,:) = frame;
%     mxs(i) = mx;
    [x,y] = function_findcenter(frame);
    
    try
        aves2(i,:,:) = frame(x-20:x+20, y-20:y+20);
    catch
        aves2(i,:,:) = nan(41,41);
    end
    
    mxs(i) = max(aves2, [], 'all');
end

bits = (2^8)-1;

figure(919)
clf
yyaxis left
% plot(exposure_times/1000, mxs/bits*100,  '-o')
plot(frames_avg, double(mxs)/bits*100, '-o')
ylabel('px val')
hold on
yyaxis right
% plot(exposure_times/1000, t*1000,  '-o')
plot(frames_avg, t*1000,  '-o')
ylabel('time (ms)')
xlabel('exposure time (ms)')
set(gca, 'XScale', 'log')
hold off

figure(929)
clf
nplts = numel(mxs);
c=0;
for i = 1:numel(mxs)
    c=c+1;
    subplot(2,nplts,c)
    imagesc(squeeze(aves(i,:,:)))
    colorbar
end
hold on
for i = 1:numel(mxs)
    c=c+1;
    subplot(2,nplts,c)
    imagesc(squeeze(aves2(i,:,:)))
    colorbar
end
hold off

%% Collect PSF
figure(1); clf
Sutter.Reference = getPosition(Sutter.obj);

position = Sutter.Reference;
moveTime=moveTo(Sutter.obj,position);

 
sz = size(Bgd);
UZ= -40:5:40;
 dataUZ = zeros([sz numel(UZ)]);
disp('Collecting PSF') 
tt = tic;
 for i=1:numel(UZ)
     position = Sutter.Reference;
     
     position(3) = position(3)+UZ(i);
     moveTime=moveTo(Sutter.obj,position);
     if i==1
            pause(1)
        else
            pause(0.1);
    end
        
         mssend(masterSocket,[pwr/1000 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        frame = thorAcquireSimple(1, 1000e3, cam);%function_BasGetFrame(Setup,3);%function_Basler_get_frames(Setup, 3 );
        frame = uint16(mean(frame,3));
        
         mssend(masterSocket,[0 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        
        frame =  max(frame-Bgd,0);
        frame = imgaussfilt(frame,2);
        dataUZ(:,:,i) =  frame;
        
        figure(1); 
        subplot(1,2,1)
        imagesc(frame);
        title(['Frame ' num2str(i)])
        subplot(1,2,2)
        imagesc(max(dataUZ,[],3));
        colorbar
        title('Max Projection')
        drawnow;
 end

        position = Sutter.Reference;
          moveTime=moveTo(Sutter.obj,position);  
disp('done')
toc(tt)

%%Calc FWHM
pxPerMu = size(frame,1) / 1000; %really rough approximate of imaging size

mxProj = max(dataUZ,[],3);
[ x,y ] =function_findcenter(mxProj );

range =2;

dimx = max((x-range),1):min((x+range),size(frame,1));
dimy =  max((y-range),1):min((y+range),size(frame,2));

thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
[a peakPlane ] = max(thisStack);
peakFrame = dataUZ(:,:,peakPlane); 

%% determine pxPerMu correctly
muUsed= 50;

disp('Determining pxPerMu')
 position = Sutter.Reference;
 position(3) = position(3)+UZ(peakPlane);
 moveTime=moveTo(Sutter.obj,position);
 
   mssend(masterSocket,[pwr/1000 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        frame = thorAcquireSimple(3, 1000e3, cam);%function_BasGetFrame(Setup,3);%function_Basler_get_frames(Setup, 3 );
        frame = uint16(mean(frame,3));
        
         mssend(masterSocket,[0 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        
        frame =  max(frame-Bgd,0);
        pos1 = imgaussfilt(frame,2);
        
        position = Sutter.Reference;
 position(3) = position(3)+UZ(peakPlane);
 position(2) = position(2)+muUsed;
 moveTime=moveTo(Sutter.obj,position);
 mssend(masterSocket,[pwr/1000 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        frame = thorAcquireSimple(3, 1000e3, cam);%function_BasGetFrame(Setup,3);%function_Basler_get_frames(Setup, 3 );
        frame = uint16(mean(frame,3));
        
         mssend(masterSocket,[0 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.01);
        end
        
        frame =  max(frame-Bgd,0);
        pos2 = imgaussfilt(frame,2);
[ x1,y1 ] =function_findcenter(pos1 );
[ x2,y2 ] =function_findcenter(pos2 );
distance = pdist([x1 y1; x2 y2]);
pxPerMu = distance/muUsed;
disp(['Determined to be ' num2str(pxPerMu) ' px per mu'])
        
xLine = double(peakFrame(x,:));
xSize = linspace(1,1000,numel(xLine))./pxPerMu; 
f1 = fit(xSize', xLine', 'gauss1');
xValue =f1.a1;
xDepth =f1.b1;
xFWHM = 2*sqrt(2*log(2))*f1.c1/sqrt(2)


ff = fit(UZ', thisStack, 'gauss1');
peakValue =ff.a1;
peakDepth =ff.b1;
peakFWHM = 2*sqrt(2*log(2))*ff.c1/sqrt(2)

figure(3);clf
subplot(1,2,1)
plot(UZ,thisStack)
hold on
plot(ff)
title(['Axial Fit FWHM ' num2str(peakFWHM) '\mum'])
subplot(1,2,2)
plot(xSize,xLine) 
hold on
plot(f1)
title(['Radial Fit FWHM ' num2str(xFWHM) '\mum'])



