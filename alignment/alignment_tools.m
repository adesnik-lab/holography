%% run this to setup SLM from scratch...

addpath(genpath('C:\Users\Holography\Documents\MATLAB\msocket\'));
rmpath(genpath('C:\Users\Holography\Documents\GitHub\SLM-Managment\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\New_SLM_Code\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\NOVOCGH_Code\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\IanTestCode\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\MSSocket_SLM\'));
%addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\SLM_MATLAB\'));
%disp establishing write protocol to master
disp('done pathing')

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.verbose = 0;
Setup.useGPU = 1;
cycleiterations = 1; % Change this number to repeat the sequence N times instead of just once
Setup.TimeToPickSequence = 0.05;    %second window to select sequence ID
Setup.SLM.timeout_ms = 1700;     %No more than 2000 ms until time out
calibID =1;                     % Select the calibration ID (z1=1 but does not exist, Z1.5=2, Z1 sutter =3);

disp('Loading Current Calibration...')
load('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\ActiveCalib.mat','CoC');
disp('Loaded.')

[ Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
[ Setup.SLM ] = Function_Start_SLM( Setup.SLM );

[Setup] = function_startBasCam(Setup);

disp('ready')
%% if you want to pair
%run this before DAQ
disp('Waiting for msocket communication From DAQ')
%then wait for a handshake
srvsock = mslisten(4211);
masterSocket = msaccept(srvsock,30);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);
%MasterIP = '128.32.177.217';
%masterSocket = msconnect(MasterIP,3002);

invar = [];

while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
end
disp('communication from Master To Holo Established');

%% blank phase on SLM
blankHolo = zeros([1920 1152]);
Function_Feed_SLM( Setup.SLM, blankHolo);
disp('sent a blank phase')

%% shoot a single holo

slmCoordsTemp = [.4 .4 -.02 1];%[0.276 .4633 .01676 1];%[0.4 0.75 0.01 1];


DEestimateTemp = DEfromSLMCoords(slmCoordsTemp); %
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimateTemp)])
[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoordsTemp );
Function_Feed_SLM( Setup.SLM, HoloTemp);
disp('sent SLM')

figure(124)
clf
imagesc(HoloTemp)
title('Hologram sent to SLM')

%% send multitarget

slmCoordsTemp = [[0.40 .40 0 1];...
                 [0.40 .60 0 1];...
                 [0.60 .60 0 1];...
                 [0.60 .40 0 1]];
                  
   
DEestimateTemp = DEfromSLMCoords(slmCoordsTemp); %
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimateTemp)])
[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,slmCoordsTemp );
Function_Feed_SLM( Setup.SLM, HoloTemp);
disp('sent SLM')

figure(124)
clf
imagesc(HoloTemp)
title('Hologram sent to SLM')

%% set laser power       
pwr =18;50;20;%;14.3; 12.5;12; %40; %13 at full; 50 at 15 divided
mssend(masterSocket,[pwr/1000 1 1]);
%% get a basler preview
function_BasPreview(Setup);
% mssend(masterSocket,[0 1 1]);

%% preview in basler multi
mssend(masterSocket,[pwr/1000 1 1]);
function_BasPreview(Setup);
mssend(masterSocket,[0 1 1]);
