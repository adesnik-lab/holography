%% Pathing and setup

makePaths()
addpath('cameras\bascam')

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
Setup.useGPU =0;
Setup.SLM.is_onek = 1;

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

disp('Ready')
%% look for objective in 1p

bas.preview()

%% connect to DAQ computer

%run this first then code on daq
disp('Waiting for msocket communication From DAQ... ')
%then wait for a handshake
srvsock = mslisten(42118);
masterSocket = msaccept(srvsock,15);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);

invar = [];

while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
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

pwr = 50;
slmCoords = [.4 .5 0 1];

disp(['Individual hologram power set to ' num2str(pwr) 'mW.'])

DEestimate = DEfromSLMCoords(slmCoords);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
Function_Feed_SLM(Setup.SLM, Holo);
mssend(masterSocket, [pwr/1000 1 1]);
bas.preview()
mssend(masterSocket, [0 1 1]);

%% Measure background noise in camera

bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

meanBgd = mean(bgd_frames, 'all');
stdBgd = std(double(bgd_frames), [], 'all');

%%Collect PSF
nframesCapture = 10;

figure(1); clf

Sutter.Reference = getPosition(Sutter.obj);

position = Sutter.Reference;
moveTo(Sutter.obj,position);

sz = size(bgd);
UZ= linspace(-70,70,21);

dataUZ = zeros([sz numel(UZ)]);

disp('Collecting PSF')

for i=1:numel(UZ)

    position = Sutter.Reference;
    position(3) = position(3)+UZ(i);
    moveTo(Sutter.obj,position);

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

    data = bas.grab(nframesCapture);

    data = mean(data, 3);
    frame = max(data-bgd, 0);
    frame = imgaussfilt(frame, 2);

    mssend(masterSocket,[0 1 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.01);
    end

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

mxProj = max(dataUZ,[],3);
[ x,y ] =function_findcenter(mxProj );

range =2;

dimx = max((x-range),1):min((x+range),size(frame,1));
dimy =  max((y-range),1):min((y+range),size(frame,2));

thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
[a peakPlane ] = max(thisStack);
peakFrame = dataUZ(:,:,peakPlane);

% %%determine pxPerMu correctly
% muUsed= 50;
% 
% disp('Determining pxPerMu')
% position = Sutter.Reference;
% position(3) = position(3)+UZ(peakPlane);
% moveTime=moveTo(Sutter.obj,position);
% 
% mssend(masterSocket,[pwr/1000 1 1]);
% invar=[];
% while ~strcmp(invar,'gotit')
%     invar = msrecv(masterSocket,0.01);
% end
% 
% data = bas.grab(nframesCapture);
% data = mean(data, 3);
% 
% mssend(masterSocket,[0 1 1]);
% invar=[];
% while ~strcmp(invar,'gotit')
%     invar = msrecv(masterSocket,0.01);
% end
% 
% frame = max(data-bgd, 0);
% pos1 = imgaussfilt(frame, 2);
% 
% position = Sutter.Reference;
% position(3) = position(3)+UZ(peakPlane);
% position(2) = position(2)+muUsed;
% moveTo(Sutter.obj,position);
% 
% mssend(masterSocket,[pwr/1000 1 1]);
% 
% invar=[];
% while ~strcmp(invar,'gotit')
%     invar = msrecv(masterSocket,0.01);
% end
% 
% data = bas.grab(nframesCapture);
% data = mean(data, 3);
% 
% mssend(masterSocket,[0 1 1]);
% invar=[];
% while ~strcmp(invar,'gotit')
%     invar = msrecv(masterSocket,0.01);
% end
% 
% frame = max(data-bgd, 0);
% pos2 = imgaussfilt(frame, 2);
% 
% [ x1,y1 ] =function_findcenter(pos1 );
% [ x2,y2 ] =function_findcenter(pos2 );
% distance = pdist([x1 y1; x2 y2]);
% pxPerMu = distance/muUsed;
% disp(['Determined to be ' num2str(pxPerMu) ' px per mu'])
% 
% %move back home
% % position = Sutter.Reference;
% % moveTime=moveTo(Sutter.obj,position);
% 
% position = Sutter.Reference;
% position(2) = position(2)+5; %add slight offset
% moveTime=moveTo(Sutter.obj,position);
% %%a
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

mimg = mean(dataUZ,3);


%%caclulations
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

