function testMultiTargetAxial
clear;close all

addpath(genpath('C:\Users\Holography\Documents\MATLAB\msocket\'));
rmpath(genpath('C:\Users\Holography\Documents\GitHub\SLM-Managment\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\New_SLM_Code\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\NOVOCGH_Code\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\Basler\'));
addpath(genpath('C:\Users\Holography\Desktop\SLM_Management\IanTestCode\'));

%disp establishing write protocol to master
disp('done pathing')


[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;

cycleiterations =1; % Change this number to repeat the sequence N times instead of just once

%Overwrite delay duration
Setup.TimeToPickSequence = 0.05;    %second window to select sequence ID
timeout = 1100;
Setup.SLM.timeout_ms = timeout;     %No more than 1100 ms until time out
calibID =1;

% load([Setup.Datapath '\07_XYZ_Calibration.mat']);
load('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\ActiveCalib.mat','CoC')

%% Generate Grid Holograms
nside = 10; %formally 13
xLoc = round(linspace(10,500,nside));
yLoc = round(linspace(10,500,nside));

testDepth =0;

SICoordinates=[];

c =0;
for i=1:numel(xLoc)
    for k = 1:numel(yLoc)
        c=c+1;
       SICoordinates(:,c)= ceil([xLoc(i) yLoc(k) testDepth]); %test Holos selected from Middle of Field zero Depth
    end
end

[SLMCoordinates] = function_SItoSLM(SICoordinates',CoC)';

AttenuationCoeffs =SLMCoordinates(4,:);
SLMCoordinates(4,:) = 1./SLMCoordinates(4,:);
SLMCoordinates(3,:) = round(SLMCoordinates(3,:),3); %Added 3/15/21 by Ian for faster compute times 

% SLMCoordinates = zeros(4,c);
% 
% SLMCoordinates(1,:) = polyvaln(COC.SI_SLM_X{calibID} ,SICoordinates');
% SLMCoordinates(2,:) = polyvaln(COC.SI_SLM_Y{calibID} ,SICoordinates');
% SLMCoordinates(3,:) = polyvaln(COC.SI_SLM_Z{calibID} ,SICoordinates');
% 
% %Add power
% disp('Calculating Diffraction Efficiency')
% SLMCoordinates(4,:) = 1./function_Power_Adjust( SLMCoordinates(1:3,:)',COC );
% AttenuationCoeffs = function_Power_Adjust( SLMCoordinates(1:3,:)', COC );


eligibleHolos = AttenuationCoeffs>0.33;
disp([num2str(sum(eligibleHolos)) ' eligible Holos']);
% tempAC =AttenuationCoeffs;
% tempAC(tempAC<0)=0;
% % singleDE = sqrt(tempAC);
% singleDE = (tempAC).^2;
singleDE=AttenuationCoeffs;
singleDE(~eligibleHolos)=[];

%%
disp('Compiling Multi Hologram')
t=tic;
Setup.verbose =0;
Setup.useGPU=1;
Setup.CGHMethod=2;


subcoordinates = SLMCoordinates(:,eligibleHolos);
[ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );

multiHologram = Hologram;

myattenuation = AttenuationCoeffs(eligibleHolos);
energy = 1./myattenuation; energy = energy/sum(energy);
multiDE = sum(energy.*myattenuation);

subcoordUnbalanced = subcoordinates;
subcoordUnbalanced(4,:)=1;
[ multiHologramUnBalanced,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordUnbalanced' );


disp(['Took  ' num2str(toc(t)) 's'])

disp('Compiling single Holograms')
hololist=[];
holoIDs= find(eligibleHolos);
for i=1:numel(holoIDs);
    t=tic;
    fprintf(['Holo ' num2str(i)]);
    subcoordinates = SLMCoordinates(:,holoIDs(i));
    [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
    hololist(:,:,i)=Hologram;
    
    %       myattenuation = AttenuationCoeffs(holoIDs(i));
    %         energy = 1./myattenuation; energy = energy/sum(energy);
    %         DE(i) = sum(energy.*myattenuation);
    %
    
    fprintf([' took ' num2str(toc(t)) 's\n']);
end
disp('done')

%% Make mSocketConnection with DAQ Computer

disp('Waiting for msocket communication')
%then wait for a handshake
srvsock = mslisten(3027);
masterSocket = msaccept(srvsock);
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


%%
try [Setup] = function_stopBasCam(Setup); catch;end
[Setup] = function_startBasCam(Setup);

% Bgdframe =  function_Basler_get_frames(Setup, 1 );
Bgdframe =  double(function_BasGetFrame(Setup,3)); %double(function_Basler_get_frames(Setup, 1 ));
Bgdframe = uint8(mean(Bgdframe,3));


[Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
[ Setup.SLM ] = Function_Start_SLM( Setup.SLM );

 Function_Feed_SLM( Setup.SLM, multiHologram);
 pause(0.02);
 disp('Ready')
%%

pwr = input('Type how much power you want to use in mW (suggest 13mW): ');
mssend(masterSocket,[pwr/1000 multiDE sum(eligibleHolos)]);
%%
% function_Basler_Preview(Setup, 5);
function_BasPreview(Setup);

mssend(masterSocket,[0/1000 multiDE sum(eligibleHolos)]);


%% Capture Multi Frame image
 Function_Feed_SLM( Setup.SLM, multiHologram);
pause(1);

mssend(masterSocket,[pwr/1000 multiDE sum(eligibleHolos)]);
invar=[];
while ~strcmp(invar,'gotit')
    invar = msrecv(masterSocket,0.1);
end
multiFrameInit =  function_BasGetFrame(Setup,1);
mssend(masterSocket,[0 1 1]);
invar=[];
while ~strcmp(invar,'gotit')
    invar = msrecv(masterSocket,0.1);
end
disp('Done')

figure(2);clf
imagesc(multiFrameInit)

 Function_Feed_SLM( Setup.SLM, multiHologramUnBalanced);
pause(1);

mssend(masterSocket,[pwr/1000 multiDE sum(eligibleHolos)]);
invar=[];
while ~strcmp(invar,'gotit')
    invar = msrecv(masterSocket,0.1);
end
multiFrameInitUnBalanced =  function_BasGetFrame(Setup,1);
mssend(masterSocket,[0 1 1]);
invar=[];
while ~strcmp(invar,'gotit')
    invar = msrecv(masterSocket,0.1);
end
disp('Done')

figure(3);clf
imagesc(multiFrameInitUnBalanced)

%% Determine the Location of Each Target

singlePwr = 12;15;%25;

clear dimx dimy
singleFrames = uint8(zeros([size(multiFrameInit) sum(eligibleHolos)]));
singleFramesCorrect = uint8(zeros([size(multiFrameInit) sum(eligibleHolos)]));

for i = 1:sum(eligibleHolos)
    disp(['Capturing From Target ' num2str(i) ' of ' num2str(sum(eligibleHolos))])
    
    Function_Feed_SLM( Setup.SLM, hololist(:,:,i));
    mssend(masterSocket,[singlePwr/1000 mean(singleDE) 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.1);
    end
    
    frame =  function_BasGetFrame(Setup,1);
    
    mssend(masterSocket,[singlePwr/1000 singleDE(i)^2 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.1);
    end
    frame2 =  function_BasGetFrame(Setup,1);

    
    mssend(masterSocket,[0 1 1]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.1);
    end
    frame =  max(frame-Bgdframe,0);
    frame = imgaussfilt(frame,2);
    
    singleFrames(:,:,i)=frame;
    
    range=15;
    
    [ x,y ] =function_findcenter(frame );
    tempx = max((x-range),1):min((x+range),size(frame,1));
    tempx2=tempx;
    tempy = max((y-range),1):min((y+range),size(frame,2));
    tempy2=tempy;
    if numel(tempx)<(range*2+1)
        tempx = padarray(tempx', (range*2+1)-numel(tempx),NaN,'post');
    end
    dimx(:,i)=tempx;
    
    if numel(tempy)<(range*2+1)
        tempy = padarray(tempy', (range*2+1)-numel(tempy),NaN,'post');
    end
    dimy(:,i)=tempy;
    
    
    smframe = frame(tempx2',tempy2');
    figure(1);
    subplot(2,3,1);
    imagesc(frame);colorbar
    subplot(2,3,2);
    imagesc(smframe);colorbar
    subplot(2,3,3);
    imagesc(hololist(:,:,i));
    
    frame2 =  max(frame2-Bgdframe,0);
    frame2 = imgaussfilt(frame2,2);
    
    singleFramesCorrect(:,:,i)=frame2;
    smframe2 = frame2(tempx2',tempy2');
    subplot(2,3,4);
    imagesc(frame2);colorbar
    subplot(2,3,5);
    imagesc(smframe2);colorbar
    
    
    drawnow;
end
disp('Done')
%%
colorLims =[10 250];
figure(2);
ax{1} = subplot(2,2,1);
imagesc(multiFrameInit);
title({[num2str(sum(eligibleHolos)) ' Targets Simultaneoulsy Displayed']; 'Power Balanced'})
colorbar
% caxis(colorLims)
ax{2} = subplot(2,2,2);
imagesc(multiFrameInitUnBalanced);
title({[num2str(sum(eligibleHolos)) ' Targets Simultaneoulsy Displayed']; 'No Power Balanced'})
colorbar
% caxis(colorLims)

ax{3} = subplot(2,2,3);
imagesc(max(singleFramesCorrect,[],3))
title({'Max image of holograms individually displayed'; 'With Balance'})
colorbar
% caxis(colorLims)
ax{4} = subplot(2,2,4);
imagesc(max(singleFrames,[],3))
title({'Max image of holograms individually displayed'; 'No power Balance'})
colorbar
linkaxes([ax{:}])
axis equal
% caxis(colorLims)

%%
imagetest(:,:,1)=double(singleBalancedAll)./double(max(singleBalancedAll(:)));
imagetest(:,:,2)=double(multiFrameInit)./double(max(multiFrameInit(:)));
figure(377);clf;image(imagetest)
axis equal

%% compute values
maxSingle = max(singleFrames,[],3);
maxSingleCorr = max(singleFramesCorrect,[],3);
ValMulti =zeros([1 size(dimy,2)]);
ValSingle =zeros([1 size(dimy,2)]);
ValSingleCorr =zeros([1 size(dimy,2)]);
for i=1:size(dimy,2);
    dx = dimx(:,i);
    dx(isnan(dx))=[];
    dy = dimy(:,i);
    dy(isnan(dy))=[];
    ValMulti(i) = mean(multiFrameInit(dx,dy),'all');
    ValMultiUB(i) = mean(multiFrameInitUnBalanced(dx,dy),'all');

    ValSingle(i) = mean(maxSingle(dx,dy),'all');
    ValSingleCorr(i) = mean(maxSingleCorr(dx,dy),'all');
end
figure(3); clf
subplot(1,2,1)
hold on
% plot(ValMulti,'o','linewidth',2,'color','k')
plot(ValMultiUB,'o','linewidth',2,'color','b')

plot(ValSingle,'o','linewidth',2,'color','r')
% plot(ValSingleCorr,'o','linewidth',2,'color',[0.5 0.5 0.5])

subplot(1,2,2);
hold on
% plot(ValMulti./max(ValMulti),'o','linewidth',2,'color','k')
% plot(ValMultiUB./max(ValMultiUB),'o','linewidth',2,'color','b')
% 
plot(ValSingle./max(ValSingle),'o','linewidth',2,'color','r')
plot(ValSingleCorr./max(ValSingleCorr),'o','linewidth',2,'color',[0.5 0.5 0.5])
%%
figure(376);clf
subplot(1,2,1)
plot(ValMultiUB./max(ValMultiUB),ValSingle./max(ValSingle),'o')
refline
subplot(1,2,2)
plot(ValMulti./max(ValMulti),ValSingleCorr./max(ValSingleCorr),'o')
refline

%%

out.singleFrames=singleFrames;
out.singleFramesCorrect=singleFramesCorrect;
out.multiFrameInit=multiFrameInit;

out.multiDE= multiDE;
out.AttenuationCoeffs = AttenuationCoeffs;
out.singleDE = singleDE;

out.SICoordinates = SICoordinates;
out.testDepth = testDepth;
out.xLoc = xLoc;
out.yLoc = yLoc;

out.ValMulti = ValMulti;
out.ValSingle = ValSingle;
out.ValSingleCorr = ValSingleCorr;
out.ValMultiUB = ValMultiUB;

save(['TestMultiPower_' date],'out')
%% Orthographic imaging
try; function_close_sutter( Sutter ); end
[ Sutter ] = function_Sutter_Start( Setup );

disp('Start with zeroing the imaging plane, show focus please')
function_BasPreview(Setup);
Sutter.Reference = getPosition(Sutter.obj);

%% Hayley z cal check
%% image multi unbalanced and single unbalanced
figure(666)
smallUZ = -30:5:30;%linspace(-20,100,50);%linspace(-60,60,50);

 Function_Feed_SLM( Setup.SLM, multiHologramUnBalanced);

 dataUZMulti = uint8(nan([size(multiFrameInit)  numel(smallUZ)])); %linspace(-30,30,15); edited 2018 02 27
dataUZSingle = uint8(zeros([size(multiFrameInit) sum(eligibleHolos) numel(smallUZ)]));%uint8(zeros([size(dimx,1) size(dimy,1) sum(eligibleHolos) numel(smallUZ)]));

 for i = 1:numel(smallUZ)
     disp(['imaging plane ' num2str(i)])
     tic
     position = Sutter.Reference;
     position(3) = position(3) + (smallUZ(i));
     moveTime=moveTo(Sutter.obj,position);
     pause(0.1);
     
     %image unbalanced multi hologram
      Function_Feed_SLM( Setup.SLM, multiHologramUnBalanced);
     mssend(masterSocket,[pwr/1000 multiDE sum(eligibleHolos)]);
     invar=[];
     while ~strcmp(invar,'gotit')
         invar = msrecv(masterSocket,0.1);
     end
    frame =  function_BasGetFrame(Setup,1);
%      frame = mean(frame,3);
     
     mssend(masterSocket,[0 1 1]);
     invar=[];
     while ~strcmp(invar,'gotit')
         invar = msrecv(masterSocket,0.1);
     end
          dataUZMulti(:,:,i) =  frame;
       figure(600);
       imagesc(frame);
       drawnow
    
       %image all the singles unbalanced
      for j = 1:sum(eligibleHolos)
        disp(['Capturing From Target ' num2str(j) ' of ' num2str(sum(eligibleHolos))])

        Function_Feed_SLM( Setup.SLM, hololist(:,:,j));
        mssend(masterSocket,[singlePwr/1000 mean(singleDE) 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.1);
        end

        frame =  function_BasGetFrame(Setup,1);

        mssend(masterSocket,[0 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.1);
        end
        frame =  max(frame-Bgdframe,0);
        frame = imgaussfilt(frame,2);
        
        tempx = dimx(:,j);
        tempx(isnan(tempx))=1;
        tempy = dimy(:,j);
        tempy(isnan(tempy))=1;
        dataUZSingle(:,:,j,i) = frame;%(tempx,tempy);
        figure(666);
        imagesc(frame(tempx,tempy)); drawnow;
      end
      toc
 end

 
 
 mssend(masterSocket,[0 1 1]);
 invar=[];
 while ~strcmp(invar,'gotit')
     invar = msrecv(masterSocket,0.1);
 end
 position = Sutter.Reference;
 moveTime=moveTo(Sutter.obj,position);
  disp('Done')
 
%% Calculate z locations multihologram 
 clear valHoloMulti peakHoloMulti valHoloSingle peakHoloSingle
  for i = 1:sum(eligibleHolos)
       tempx = dimx(:,i);
     tempx(isnan(tempx))=[];
     tempy = dimy(:,i);
     tempy(isnan(tempy))=[];
     
     holoDat = squeeze(mean(mean(dataUZMulti(tempx,tempy,:))));
     valHoloMulti(:,i) = holoDat;
     singleDat = squeeze(mean(mean(dataUZSingle(:,:,i,:))));
     valHoloSingle(:,i) = singleDat;
     try
        ff = fit(smallUZ', squeeze(holoDat), 'gauss1');
        peakHoloMulti(i) = ff.b1;
        ff = fit(smallUZ', squeeze(singleDat), 'gauss1');
        peakHoloSingle(i) = ff.b1;
     catch
         disp(['Error on holo ' num2str(i)])
              peakHoloMulti(i) = 0;
              peakHoloSingle(i) = 0;
     end
  end
  %%
  peakHoloMultiBackup =peakHoloMulti;
  peakHoloSingleBackup=peakHoloSingle;
  %%
  peakHoloSingle = peakHoloSingleBackup;
  peakHoloMulti = peakHoloMultiBackup;
  
  peakHoloSingle(max(valHoloSingle)<.0005)=nan;
  peakHoloMulti(max(valHoloMulti)<.0005)=nan;
  
  disp([num2str(sum(max(valHoloSingle)<5 | max(valHoloMulti)<5)) ' points excluded because val too low'])
  %% quick analysis
  figure(35);clf
  plot(peakHoloMultiBackup,peakHoloSingleBackup,'o')
  ylabel('Single Target Z Loc')
  xlabel('Multi Target Z Loc')
  refline
  axis equal
  
%%
%First image the multihologram
Sutter.Reference = getPosition(Sutter.obj);
smallUZ = -20:5:80;%linspace(-20,100,50);%linspace(-60,60,50);

 Function_Feed_SLM( Setup.SLM, multiHologram);

 dataUZMulti = uint8(nan([size(multiFrameInit)  numel(smallUZ)])); %linspace(-30,30,15); edited 2018 02 27
 for i = 1:numel(smallUZ)
     position = Sutter.Reference;
     position(3) = position(3) + (smallUZ(i));
     moveTime=moveTo(Sutter.obj,position);
     pause(0.1);
     
     mssend(masterSocket,[pwr/1000 multiDE sum(eligibleHolos)]);
     invar=[];
     while ~strcmp(invar,'gotit')
         invar = msrecv(masterSocket,0.1);
     end
    frame =  function_BasGetFrame(Setup,1);
%      frame = mean(frame,3);
     
     mssend(masterSocket,[0 1 1]);
     invar=[];
     while ~strcmp(invar,'gotit')
         invar = msrecv(masterSocket,0.1);
     end
     
          dataUZMulti(:,:,i) =  frame;

          figure(1);
          subplot(1,2,1);
        imagesc(frame); title('Current Frame')
        subplot(1,2,2);
        title('Mean Image')
        imagesc(nanmean(dataUZ,3));
        drawnow
 end

 
 
 mssend(masterSocket,[0 1 1]);
 invar=[];
 while ~strcmp(invar,'gotit')
     invar = msrecv(masterSocket,0.1);
 end
 position = Sutter.Reference;
 moveTime=moveTo(Sutter.obj,position);
  disp('Done')
 
  
 %%
 dummy = input('Now Turn on ScanImage. (press enter when ready)');


dataSI = uint8(nan([size(multiFrameInit)  numel(smallUZ)])); %linspace(-30,30,15); edited 2018 02 27
for i = 1:numel(smallUZ)
    position = Sutter.Reference;
    position(3) = position(3) + (smallUZ(i));
    moveTime=moveTo(Sutter.obj,position);
    pause(0.1);
    
    frame =  function_BasGetFrame(Setup,3);
        frame = mean(frame,3);
        
    dataSI(:,:,i) =  frame;
    
    figure(1);
    subplot(1,2,1);
    imagesc(frame);
    subplot(1,2,2);
    imagesc(nanmean(dataSI,3));
    drawnow
    
end

 position = Sutter.Reference;
 moveTime=moveTo(Sutter.obj,position);
 disp('Done')

 
 %% Calculate PSFs
 disp('Fitting Gaussians...')
 clear valHolo valSI
 for i = 1:sum(eligibleHolos)
     fprintf([num2str(i) ' ']);
     tempx = dimx(:,i);
     tempx(isnan(tempx))=[];
     tempy = dimy(:,i);
     tempy(isnan(tempy))=[];
     
     holoDat = mean(mean(dataUZ(tempx,tempy,:)));
     valHolo(:,i) = holoDat;
     try
     ff = fit(smallUZ', squeeze(holoDat), 'gauss1');
     peakHolo(i) =ff.b1;
     catch
              peakHolo(i) =0;
     end

     
     SIDat =  mean(mean(dataSI(tempx,tempy,:)));
     valSI(:,i) = SIDat;
     try
     ff = fit(smallUZ', squeeze(SIDat), 'gauss1');
     peakSI(i) = ff.b1;
     catch
              peakSI(i) = 0;
     end

 end
 disp('done')
 
%  [a peakHolo]=max(valHolo);
%  [a peakSI]=max(valSI);
 

%  ff = fit(smallUZ', valHolo(:,1), 'gauss1');
 
maxSI = max(valSI);
maxHolo = max(valHolo);
%%
 figure(3);clf
 subplot(1,3,1)
 plot(valHolo)
 title('Hologram PSFs')
 
 subplot(1,3,2)
 plot(valSI)
 title('Scanimage PSFs')
 
 subplot(1,3,3)
 plot(peakHolo,peakSI,'o')
 xlabel('Peak Holo')
 ylabel('PeakSI')
 refline(1)
 
 %%
 figure(4);clf
 subplot(1,2,1)
 coords = SICoordinates(1:2,eligibleHolos);
 
 eligible = maxHolo>10 & maxHolo<200;
 
 coords = coords(:,eligible);
 pSI = peakSI(eligible);
 pHolo = peakHolo(eligible);
 
 scatter3(coords(1,:),coords(2,:),pSI,'filled')
 hold on
  scatter3(coords(1,:),coords(2,:),pHolo,'filled','r')

  title(['Optotune ' num2str(testDepth)]);
  subplot(1,2,2)
   scatter3(coords(1,:),coords(2,:),pSI-pHolo)
   zlabel('Error \mum')

   
   
  
  %%
  
  out.peakSI = peakSI;
  out.peakHolo = peakHolo;
  out.coords = coords;
  out.valSI =valSI;
  out.valHolo=valHolo;
%   out.multiFrameInit=multiFrameInit;
%   out.singleFrames = singleFrames;
  out.maxSI=maxSI;
  out.maxHolo=maxHolo;
  
 save(['TestMultiAxial_' date],'out')

  %%
  try; function_close_sutter( Sutter ); end

  