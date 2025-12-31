%% Pathing and setup

makePaths()

clear
close all
clc

tBegin = tic;

disp('Setting up stuff...');

Setup = function_loadparameters3();

if Setup.useGPU
    disp('Getting gpu...'); %this can sometimes take a while at initialization
    g= gpuDevice;
end

slm = MeadowlarkSLM();
slm.start()

disp('SLM Ready!')
%% Sutter

Setup.Sutterport ='COM3';
try; function_close_sutter( Sutter ); end
[ Sutter ] = function_Sutter_Start( Setup );

disp('Sutter Ready!')

%% Basler

bas = bascam();
bas.start()

disp('Basler Ready!')

%% Basler preview

bas.preview()

%% optional:  connect to DAQ computer

%run this first then code on daq
disp('Waiting for msocket communication From DAQ... ')
%then wait for a handshake
srvsock = mslisten(42136);
masterSocket = msaccept(srvsock,15);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);

invar = [];

while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
end
fprintf('done.\r')

%% blank phase on SLM

slm.blank()

%% shoot single holo, no power control

slmCoordsTemp = [0.55 0.55 0 1];

[ HoloTemp,~,~ ] = function_Make_3D_SHOT_Holos(Setup, slmCoordsTemp);
DEestimateTemp = DEfromSLMCoords(slmCoordsTemp);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimateTemp)])
slm.feed(HoloTemp)
disp('sent SLM')

figure(124)
clf
imagesc(HoloTemp)
colorbar
title('Hologram sent to SLM')

%% shoot multi holo, no power control

slmCoordsTemp = [[0.30 .40 .10 1];...
                 [0.45 .66 .2 1];...
                 [0.78 .77 -.03 1];...
                 [0.60 .40 0 1]];

[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos(Setup, slmCoordsTemp);
% Function_Feed_SLM(Setup.SLM, HoloTemp);
disp('sent SLM')

figure(121)
clf
imagesc(HoloTemp)
title('Hologram sent to SLM')

%% shoot hologram with power control
% must be connected to daq computer

% input power in mW
pwr = 250;

% slmCoordsTemp = [0.55 0.55 0.015 1];
slmCoordsTemp = [[0.13 .15 0.02 1];...
                 [0.05 .80 0.02 1];...
                 [1.0 .4 0.02 1];...
                 [0.9 .9 0.02 1];...
                 [0.60 .60 0.02 1]];

% slmCoordsTemp = [[0.20 .17 0.02 1];... % top right
%                  [0.20 .80 0.02 1];... % bottom right
%                  [0.79 .20 0.02 1];... % top left
%                  [0.81 .83 0.02 1];...
%                  [0.60 .60 0.02 1]];

[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos(Setup, slmCoordsTemp);
slm.feed(HoloTemp)
disp('sent SLM')
figure(124)
imagesc(HoloTemp)

mssend(masterSocket, [pwr/1000 1 1]);
bas.preview()
% figure
% frame = bas.grab(3);
% imagesc(mean(frame, 3), [0, 24])
mssend(masterSocket, [0 1 1]);
