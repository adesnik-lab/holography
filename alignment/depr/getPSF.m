function getPSF

%% Pathing
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

disp('starting basler...')
try function_stopBasCam(Setup); end
[Setup] = function_startBasCam(Setup);
disp('Ready')
%% look for objective in 1p or you know... have it already set up

function_BasPreview(Setup);

%% Make mSocketConnections with DAQ and SI Computers

disp('Waiting for msocket communication From DAQ')
%then wait for a handshake
srvsock = mslisten(42163);
masterSocket = msaccept(srvsock,10);
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

pwr = 25;%;14.3; 12.5;12; %40; %13 at full; 50 at 15 divided
disp(['individual hologram power set to ' num2str(pwr) 'mW']);

%%
disp('Find the spot and check if this is the right amount of power')
slmCoordsTemp =[.45 .55 -.05 1];%[0.276 .4633 .01676 1];%[0.4 0.75 0.01 1];
DEestimateTemp = DEfromSLMCoords(slmCoordsTemp); %
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimateTemp)])

[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoordsTemp );

blankHolo = zeros([1920 1152]);
%Function_Feed_SLM( Setup.SLM, blankHolo);

Function_Feed_SLM( Setup.SLM, HoloTemp);

mssend(masterSocket,[pwr/1000 1 1]);

function_BasPreview(Setup);
mssend(masterSocket,[0 1 1]);

%% Collect background frames for signal to noise testing
disp('Collecting Background Frames');

nBackgroundFrames = 10;
nFramesCapture = 10;

Bgdframe = function_BasGetFrame(Setup,nBackgroundFrames);% function_Basler_get_frames(Setup, nBackgroundFrames );
Bgd = uint16(mean(Bgdframe,3));
meanBgd = mean(single(Bgdframe(:)));
stdBgd =  std(single(Bgdframe(:)));

threshHold = meanBgd+3*stdBgd;

fprintf(['3\x03c3 above mean threshold ' num2str(threshHold,4) '\n'])

%%Collect PSF
figure(1); clf
Sutter.Reference = getPosition(Sutter.obj);

position = Sutter.Reference;
moveTime=moveTo(Sutter.obj,position);

 
sz = size(Bgd);
UZ= linspace(-40,40,13);%-40:5:40;
 dataUZ = zeros([sz numel(UZ)]);
disp('Collecting PSF') 
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
        frame = function_BasGetFrame(Setup,nFramesCapture);%function_Basler_get_frames(Setup, 3 );
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

%%Calc FWHM
pxPerMu = size(frame,1) / 2000; %really rough approximate of imaging size

mxProj = max(dataUZ,[],3);
[ x,y ] =function_findcenter(mxProj );

range =2;

dimx = max((x-range),1):min((x+range),size(frame,1));
dimy =  max((y-range),1):min((y+range),size(frame,2));

thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
[a peakPlane ] = max(thisStack);
peakFrame = dataUZ(:,:,peakPlane); 

%%determine pxPerMu correctly
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
        frame = function_BasGetFrame(Setup,nFramesCapture);%function_Basler_get_frames(Setup, 3 );
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
        frame = function_BasGetFrame(Setup,10);%function_Basler_get_frames(Setup, 3 );
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

%move back home
% position = Sutter.Reference;
% moveTime=moveTo(Sutter.obj,position);

position = Sutter.Reference;
position(2) = position(2)+5; %add slight offset
moveTime=moveTo(Sutter.obj,position);
        %%a
xLine = double(peakFrame(x,:));
xSize = linspace(1,1000,numel(xLine))./pxPerMu; 
f1 = fit(xSize', xLine', 'gauss1');
xValue =f1.a1;
xDepth =f1.b1;
xFWHM = 2*sqrt(2*log(2))*f1.c1/sqrt(2);


ff = fit(UZ', thisStack, 'gauss1');
peakValue =ff.a1;
peakDepth =ff.b1;
peakFWHM = 2*sqrt(2*log(2))*ff.c1/sqrt(2);

results.UZ = UZ;
results.thisStack = thisStack;
results.ff = ff;
results.peakFWHM = peakFWHM;
results.xSize = xSize;
results.xLine = xLine;
results.f1 = f1;
results.xFWHM = xFWHM;
mimg = mean(dataUZ,3);
results.mimg = mimg;
results.peakFrame = peakFrame;
results.peakVal =max(peakFrame(:)); 


figure(3);clf
subplot(2,2,1)
plot(UZ,thisStack)
hold on
plot(ff)
title(['Axial Fit FWHM ' num2str(peakFWHM) '\mum'])
ylabel('Intensity')
xlabel('Z-position (\mum)')
legend('Measured', 'Fit')

subplot(2,2,2)
plot(xSize,xLine) 
hold on
plot(f1)
xlim([f1.b1-50 f1.b1+50])
title(['Radial Fit FWHM ' num2str(xFWHM) '\mum'])
ylabel('Intensity')
xlabel('X-position (\mum)')
legend('Measured', 'Fit')



subplot(2,2,3);
axis square
imagesc(mimg)
% imagesc(mxProj)
xlim([y-30 y+30])
ylim([x-30 x+30])
title('Mean Holo')
colorbar

subplot(2,2,4);
axis square
imagesc(peakFrame)
% imagesc(mxProj)
xlim([y-30 y+30])
ylim([x-30 x+30])
title('Peak Frame')
colorbar

max(peakFrame(:))
%%

dataOut.p276g = results