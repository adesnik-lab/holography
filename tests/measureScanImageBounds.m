%% Pathing and setup


clear
close all
clc

tBegin = tic;

disp('Setting up stuff...');

makePaths2()

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

sutter = sutterController();

bas = bascam();
bas.start()

disp('Ready')

%% connect to DAQ computer

%run this first then code on daq
fprintf('Waiting for msocket communication From DAQ... ')
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

%% align the ScanImage into cam and find Z plane
bas.preview()

%% measure scanimage scan boundary

data = bas.grab(10);
data = mean(data,3);

scale = 10;
data= imresize(data, 1/scale);

figure(1)
clf
imagesc(data)
colorbar
%%
thresh = 4;

figure(2)
clf
subplot(1,2,1)
imagesc(data)
colorbar

subplot(1,2,2)
im = imgaussfilt(data,2);
im = im>thresh;
imagesc(im)
colorbar

[r,c] = find(im);

SIboundary = boundary(r,c);
SIxboundary = c(SIboundary);
SIyboundary = r(SIboundary);
hold on
p=plot(SIxboundary,SIyboundary,'r');


%% find holo on plane

pwr = 10;
slmCoords = [0.6 0.6 -0.035 1];

disp(['Individual hologram power set to ' num2str(pwr) 'mW.'])

DEestimate = DEfromSLMCoords(slmCoords);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
Function_Feed_SLM(Setup.SLM, Holo);
mssend(masterSocket, [pwr/1000 1 1]);


f = figure('Name','Basler Preview', 'NumberTitle','off');
while isgraphics(f)
    frame = bas.grab(1);
%     frame= imresize(frame, 1/scale);
    imagesc(frame)
%     caxis([0 2^bas.bit_depth-1])
    colorbar
    hold on
    p=plot(SIxboundary*scale,SIyboundary*scale,'r', 'LineWidth',2);
    drawnow
end


mssend(masterSocket, [0 1 1]);







