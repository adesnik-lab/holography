timeout = 1700;

addpath(genpath('C:\Users\Holography\Desktop\holography'))
addpath(genpath('C:\Users\Holography\Desktop\meadowlark'))

Setup = function_loadparameters3();
Setup.CGHMethod = 2; % now defaults to GSS
Setup.verbose = 0;
% Setup.useGPU = 1; % now defaults to GPU

cycleiterations = 1; % Change this number to repeat the sequence N times instead of just once

%Overwrite delay duration
Setup.TimeToPickSequence = 0.05;    %second window to select sequence ID
Setup.SLM.timeout_ms = timeout;     %No more than 2000 ms until time out
% calibID = 1;                        % Select the calibration ID (z1=1 but does not exist, Z1.5=2, Z1 sutter =3);
Setup.calib = 'C:\Users\Holography\Desktop\calibs\ActiveCalib.mat';

disp('Loading current calibration...')
load(Setup.calib,'CoC');
disp(['Successfully loaded CoC from: ', Setup.calib])
%% now start msocket communication
 
disp('Waiting for msocket communication')

srvsock = mslisten(42118); %3027
masterSocket = msaccept(srvsock,6);
msclose(srvsock)
sendVar = 'A';
mssend(masterSocket, sendVar);

invar = [];
while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
end

disp('communication from Master To Holo Established')

x = 1;     
HRin = []; 
while x>0
    HRin = msrecv(masterSocket,.5);

    if ~isempty(HRin);
        disp('new File Detected - running HoloRequest')
        holoRequest = HRin;
        x=0;
    end
end

%%Load SI Coordinates

if ~isfield(holoRequest,'ignoreROIdata')  %if we're doing things the old way
    try
        load([Setup.Holorequestpath 'ROIData.mat']);
    catch
        disp('No ROIData file')
        return
    end
    LN = numel(ROIdata.rois);
    SICoordinates = zeros(3,LN);
    for i = 1:LN
        u = mean(ROIdata.rois(i).vertices);
        u(1)=u(1)+holoRequest.xoffset;
        u(2)=u(2)+holoRequest.yoffset;
        
        SICoordinates(1:2,i) = u;
        SICoordinates(3,i) = ROIdata.rois(i).OptotuneDepth;
    end
    SLMCoordinates = zeros(4,LN);
    
else  %if I'm doing a custom sequence
  LN = size(holoRequest.targets,1);
    SLMCoordinates = zeros(4,LN);  
    SICoordinates = holoRequest.targets;
    SICoordinates(:,1)=SICoordinates(:,1)+holoRequest.xoffset;
    SICoordinates(:,2)=SICoordinates(:,2)+holoRequest.yoffset;    
    SICoordinates=SICoordinates';
end
%%quickly compute DEs and return them over msocket

if isfield(holoRequest,'roiWeights')
    weightsToUse = holoRequest.roiWeights;
    weightsToUse(isnan(weightsToUse))=1;
    disp('Weighting Holograms based on roiWeights')
else
    weightsToUse = ones([1 LN]);
    disp('NO weights detected using flat weight')
end

[AC, DE_list] = computeDEfromList(SICoordinates, holoRequest.rois, weightsToUse);

mssend(masterSocket,DE_list);  
disp('Sent DE to master');

%%Compute SLM Coordinates
DEfloor = 0.05;

[SLMCoordinates] = function_SItoSLM(SICoordinates',CoC)';
AttenuationCoeffs =SLMCoordinates(4,:);
lowDE = AttenuationCoeffs<DEfloor;
AttenuationCoeffs(lowDE)=DEfloor;
SLMCoordinates(4,lowDE)=DEfloor;
disp([num2str(sum(lowDE)) ' Target(s) below Diffraction Efficiency floor (' num2str(DEfloor) ').']);

SLMCoordinates(4,:) = 1./SLMCoordinates(4,:);
SLMCoordinates(3,:) = round(SLMCoordinates(3,:),3); %Added 2/24/21 by Ian for faster compute times 

%%x
close
f = figure('units','normalized','outerposition',[0.125 0.5 0.75 0.5]);

subplot(1,3,1)
scatter3(SICoordinates(1,:),SICoordinates(2,:),SICoordinates(3,:),[],SLMCoordinates(4,:),'filled'); 
% colorbar;
xlabel('X, SI coordinates');ylabel('Y, SI coordinates'); zlabel('Z, SI coordinates'); title('Intensity Correction coefficients');
subplot(1,3,2)
scatter3(SLMCoordinates(1,:),SLMCoordinates(2,:),SLMCoordinates(3,:),[],SLMCoordinates(4,:),'filled'); 
% colorbar;
xlabel('X, SLM coordinates');ylabel('Y, SLM coordinates'); zlabel('Z, SLM coordinates'); title('Intensity Correction coefficients');
subplot(1,3,3)
hist(AttenuationCoeffs,20);
ylabel('Count')
xlabel('Single Target Diffraction Efficiency')
title('Histogram of Diffraction Efficiencies')
pause(1);



%%Sort Holograms

% holoRequest.rois;

numTargets = cellfun(@(x) numel(x),holoRequest.rois);

numSolo = sum(numTargets==1);
solos = find(numTargets==1);
numMid  = sum(numTargets>1 & numTargets < 25);
mids = find(numTargets>1 & numTargets < 25);
numLarge = sum(numTargets>=25 & numTargets < 50);
larges = find(numTargets>=25 & numTargets < 50);
numExtraLarge = sum( numTargets >= 50);
extraLarges =find( numTargets >= 50);


%%Compile Holograms
% optomized for our hologram computer that has a certain amount of gpu space
if size(weightsToUse,2)==1
    weightsToUse=weightsToUse';
end

SLMCoordinates(4,:)=SLMCoordinates(4,:).*weightsToUse; %I don't think this will do anything for single target holos but i'm not sure.

hololist=[];
%solos
if numSolo ==0
    disp('No Single Target Holos')
elseif numSolo<40
    disp('Less than 40 Single Holo Targets')
    for i=1:numSolo
        j=solos(i);
        disp(['Now compiling hologram ' int2str(i) ' of ' int2str(numel(solos))])
        ROIselection = holoRequest.rois{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; 
        energy = energy/sum(energy);
        DE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*DE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        hololist(:,:,j) = Hologram;
    end
else 
    disp('More than 40 Single Holo Targets')
    
    clear tempDE temphololist
    ROIs =  holoRequest.rois([solos]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=4
        delete(p);
        parpool(4);
    end
    %7/24 hayley switched to not parfor
    parfor j=1:numSolo
        disp(['Now compiling hologram ' int2str(j) ' of ' int2str(numel(solos))])
        ROIselection =ROIs{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        tempDE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*tempDE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        temphololist(:,:,j) = Hologram;
    end
    for i = 1:numSolo
        j=solos(i);
        DE(j)=tempDE(i);
        hololist(:,:,j)=temphololist(:,:,i);
    end
end 


%Midsize holograms
if numMid ==0
    disp('No Midsized Holos')
else 
    disp('Computing Midsized Holos')
    
    clear tempDE temphololist
    ROIs =  holoRequest.rois([mids]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=4
        delete(p);
        parpool(4);
    end
    parfor j=1:numMid
        disp(['Now compiling hologram ' int2str(j) ' of ' int2str(numel(mids))])
        ROIselection =ROIs{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        tempDE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*tempDE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        temphololist(:,:,j) = Hologram;
    end
    for i = 1:numMid
        j=mids(i);
        DE(j)=tempDE(i);
        hololist(:,:,j)=temphololist(:,:,i);
    end
end 


%Large holograms
if numLarge ==0
    disp('No Large Holos')
else 
    disp('Computing Large Holos')
    
    clear tempDE temphololist
    ROIs =  holoRequest.rois([larges]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=2;
        delete(p);
        parpool(2);
    end
    parfor j=1:numLarge
        disp(['Now compiling hologram ' int2str(j) ' of ' int2str(numel(larges))])
        ROIselection =ROIs{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        tempDE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*tempDE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        temphololist(:,:,j) = Hologram;
    end
    for i = 1:numLarge
        j=larges(i);
        DE(j)=tempDE(i);
        hololist(:,:,j)=temphololist(:,:,i);
    end   
end

%Extra Large holograms
if numExtraLarge ==0
    disp('No Extra Large Holos')
else 
    disp('Computing Extra Large Holos')
        
    for i=1:numExtraLarge
        j=extraLarges(i);
        disp(['Now compiling hologram ' int2str(i) ' of ' int2str(numel(extraLarges))])
        ROIselection = holoRequest.rois{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        DE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*DE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        hololist(:,:,j) = Hologram;
    end
end 

disp('Done')
%%Shoot



%locations=FrankenScopeRigFile();
%save('Y:\holography\FrankenRig\HoloRequest-DAQ\HoloRequest.mat','DE_list','-append')
%totally remove sequences as a thing that exists basically
sequences = uint8(hololist); %shouldn't change anything added 9/14/21

flushMSocket(masterSocket)

[Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
Setup.SLM.wait_For_Trigger = 1;
Setup.SLM.timeout_ms = timeout;
[ Setup.SLM ] = Function_Start_SLM( Setup.SLM );

%sendVar='B';
%mssend(masterSocket,sendVar)
orderBackup=[]; %Sequence list is archived in case the daq errors. normally disposed of after exp. 1/19/21
c=1;
while true
    orderBackup{c} = ShootSequencesMsocket(Setup,sequences,masterSocket);
    c=c+1;
end


