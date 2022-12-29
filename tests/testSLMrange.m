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
Setup.useGPU =1;
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
addpath('cameras\bascam')
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

pwr = 200;
slmCoords = [.6 .4 0 1];

disp(['Individual hologram power set to ' num2str(pwr) 'mW.'])

DEestimate = DEfromSLMCoords(slmCoords);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
Function_Feed_SLM(Setup.SLM, Holo);
mssend(masterSocket, [pwr/1000 1 1]);
bas.preview()
mssend(masterSocket, [0 1 1]);

%% TEST CODE

z_val = 0;

x_test = 0:0.1:1;
y_test = 0:0.1:1;

x_hold = 0.4;
y_hold = 0.4;

nframesCapture = 10;


%%

bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

meanBgd = mean(bgd_frames, 'all');
stdBgd = std(double(bgd_frames), [], 'all');

sz = size(bgd);

xs = [];
xUZ = zeros([numel(x_test) sz]);

for i=1:numel(x_test)
    xtest = x_test(i);
    slmCoords = [xtest y_hold z_val 1];
    [Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
    Function_Feed_SLM(Setup.SLM, Holo);
    mssend(masterSocket, [pwr/1000 1 1]);

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

    xUZ(i,:,:) = frame;
    
    [ x1,y1 ] =function_findcenter(frame );

    xs = [xs x1];
end

%%

figure(1)
clf

for i=1:11
    imagesc(squeeze(xUZ(i,:,:)))
    pause
end

%%



ys = [];













