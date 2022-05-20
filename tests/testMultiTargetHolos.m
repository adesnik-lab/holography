function testMultiTargetHolos
try [Setup] = function_stopBasCam(Setup); catch;end

clear;close all
timeout = 1100;

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
Setup.SLM.timeout_ms = timeout;     %No more than 2000 ms until time out
calibID =1;

% load([Setup.Datapath '\07_XYZ_Calibration.mat']);
load('C:\Users\Holography\Desktop\SLM_Management\Calib_Data\ActiveCalib.mat','CoC')


[Setup] = function_startBasCam(Setup);
disp('Ready')
%% Create Dummy Targets

disp('Generating Random Targets')
throwout = rand; %get rid of weird matlab rand initialization error

nTestHolos = 100;
testDepth = 0;

TN = 10000; %total Distractor Targets
depths = [0 35 55];

SICoordinates= zeros(3,TN);

for i=1:nTestHolos
    SICoordinates(:,i)= ceil([rand*255+128 rand*255+128 testDepth]); %test Holos selected from Middle of Field zero Depth
end

for i=nTestHolos:TN
    SICoordinates(:,i)= ceil([rand*511 rand*511 depths(ceil(rand*numel(depths)))]);
end

% SLMCoordinates = zeros(4,TN);
%
% SLMCoordinates(1,:) = polyvaln(COC.SI_SLM_X{calibID} ,SICoordinates');
% SLMCoordinates(2,:) = polyvaln(COC.SI_SLM_Y{calibID} ,SICoordinates');
% SLMCoordinates(3,:) = polyvaln(COC.SI_SLM_Z{calibID} ,SICoordinates');

%Add power
disp('Calculating Diffraction Efficiency')
% SLMCoordinates(4,:) = 1./function_Power_Adjust( SLMCoordinates(1:3,:)',COC );
% AttenuationCoeffs = function_Power_Adjust( SLMCoordinates(1:3,:)', COC );

[SLMCoordinates] = function_SItoSLM(SICoordinates',CoC)';
AttenuationCoeffs =SLMCoordinates(4,:);
SLMCoordinates(4,:) = 1./SLMCoordinates(4,:);
SLMCoordinates(3,:) = round(SLMCoordinates(3,:),3); %Added 2/24/21 by Ian for faster compute times 


eligibleHolos = AttenuationCoeffs>0.33;
disp([num2str(sum(eligibleHolos)) ' Targets with DE > 0.33']);
disp([num2str(sum(eligibleHolos(1:nTestHolos))) ' of ' num2str(nTestHolos) ' Test targets eligible'] );

%% Make List Holograms
tMake=tic;

loadOldHoloList=1;
if loadOldHoloList
    [f p] =uigetfile;
    load(fullfile(p,f),'out');
    
    hololist= out.holoList;
    nTargetsList = out.nTargetsList;
    holoIDList= out.holoIDList;
    DElist= out.DElist;
    containsTest= out.containsTest;
    DEeach= out.DEeach;
    SICoordinates = out.SICoordinates;
    SLMCoordinates = out.SLMCoordinates;
    
else
    disp('Making a List of Holos')
    testHolos = find(eligibleHolos(1:nTestHolos));
    distractHolos = find(eligibleHolos);
    distractHolos(1:numel(testHolos))=[];
    numDistractTotal = numel(distractHolos);
    
    nTestHoloToUse = 10; %10
    nDistractHolos = [1 2 3 5 10 25 50];%[5 10 25];
    distractIter = 5; %5
    
    nTotalHolograms = nTestHoloToUse *numel(nDistractHolos) * distractIter *2 +nTestHoloToUse ;
    
    hololist = zeros(Setup.Nx,Setup.Ny, nTotalHolograms,'uint8');
    
    clear nTargetsList holoIDList  DElist containsTest DEeach
    
    idx = 0;
    for i =1:nTestHoloToUse
        idx=idx+1;
        nTargetsList(idx) = 1;
        thisTest = testHolos(i);
        holoIDList{idx}= thisTest;
        DElist{idx}= AttenuationCoeffs(thisTest);
        myattenuation = AttenuationCoeffs(thisTest);
        energy = 1./myattenuation; energy = energy/sum(energy);
        DEeach(idx) = sum(energy.*myattenuation);
        %
        containsTest(idx) = i;
        for k=1:numel(nDistractHolos)
            for L = 1:distractIter;
                %without Test
                idx=idx+1;
                numDistractors = nDistractHolos(k);
                nTargetsList(idx) = numDistractors;
                theseDistract = distractHolos(randperm(numDistractTotal,numDistractors));
                holoIDList{idx} = theseDistract';
                DElist{idx} = AttenuationCoeffs(theseDistract');
                myattenuation = AttenuationCoeffs(theseDistract');
                energy = 1./myattenuation; energy = energy/sum(energy);
                DEeach(idx) = sum(energy.*myattenuation);
                containsTest(idx)=0;
                
                %with test
                idx=idx+1;
                
                nTargetsList(idx) = numDistractors+1;
                holoIDList{idx} = [theseDistract thisTest];
                DElist{idx} = AttenuationCoeffs([theseDistract thisTest]);
                myattenuation = AttenuationCoeffs([theseDistract thisTest]);
                energy = 1./myattenuation; energy = energy/sum(energy);
                DEeach(idx) = sum(energy.*myattenuation);
                containsTest(idx)=i;
            end
        end
    end
    
    disp([num2str(nTotalHolograms) ' Holograms selected'])
    
    %%Compile Holograms
    Setup.verbose =0;
    Setup.useGPU=1;
    Setup.CGHMethod=2;
    
      
    tAllCompile = tic;
    
    numTargets = cellfun(@(x) numel(x),holoIDList);

numSolo = sum(numTargets==1);
solos = find(numTargets==1);
numMid  = sum(numTargets>1 & numTargets < 25);
mids = find(numTargets>1 & numTargets < 25);
numLarge = sum(numTargets>=25 & numTargets < 50);
larges = find(numTargets>=25 & numTargets < 50);
numExtraLarge = sum( numTargets >= 50);
extraLarges =find( numTargets >= 50);
        hololist=[];
%solos
if numSolo ==0
    disp('No Single Target Holos')
elseif numSolo<40
    disp('Less than 40 Single Holo Targets')
    for i=1:numSolo
        j=solos(i);
        disp(['Now compiling hologram ' int2str(i) ' of ' int2str(numel(solos))])
        ROIselection = holoIDList{j};
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
    ROIs = holoIDList([solos]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=4
        delete(p);
        parpool(4);
    end
    
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
    ROIs =  holoIDList([mids]);
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
    ROIs = holoIDList([larges]);
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
        ROIselection = holoIDList{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        DE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*DE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        hololist(:,:,j) = Hologram;
    end
end 

disp(['Done. took ' num2str(toc(tAllCompile)) 's'])
%           try
%         p = gcp('nocreate');
%         if ~isfield(p,'NumWorkers')
%             parpool(4);
%         elseif p.NumWorkers ~=4
%             delete(gcp('nocreate'));
%             parpool(4);
%         end
%         catch
%              delete(gcp('nocreate'));
%             parpool(4);
%         end
%     
%         parfor j = 1:nTotalHolograms
%             %     for j = 1:nTotalHolograms
%             
%             t=tic;
%             disp(['Now compiling hologram ' int2str(j) ' of ' int2str(nTotalHolograms)])
%             ROIselection = holoIDList{j};
%             myattenuation = AttenuationCoeffs(ROIselection);
%             energy = 1./myattenuation; energy = energy/sum(energy);
%             DE(j) = sum(energy.*myattenuation);
%             %disp(['Hologram has ' num2str(numel(ROIselection)) ' Target(s)']);
%             disp(['DE of hologram with ' num2str(numel(ROIselection)) ' Target(s): ' int2str(100*DE(j)) '%']);
%             subcoordinates = SLMCoordinates(:,ROIselection);
%             [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
%             hololist(:,:,j) = Hologram;
%             fprintf(['Took ' num2str(toc(t)) 's\n']);
%         end
%     
    
    disp('Saving Hololist')
    out.holoList = hololist;
    out.nTargetsList =nTargetsList;
    out.holoIDList =holoIDList;
    out.DElist = DElist;
    out.containsTest = containsTest;
    out.DEeach = DEeach;
    out.SICoordinates = SICoordinates;
    out.SLMCoordinates = SLMCoordinates;
    
    save(['tempHololist_' date],'out','-v7.3');
end

disp(['Done. Took ' num2str(toc(tMake)) 's']);
%% Recompute attenuations (not normally nescessary
[tempSLM] = function_SItoSLM(SICoordinates',CoC)';
AttenuationCoeffs =tempSLM(4,:);
for i=1:numel(holoIDList)
    holos = holoIDList{i} ;
    DElist{i} = AttenuationCoeffs(holos');
    
    
    myattenuation = DElist{i};
    energy = 1./myattenuation;
    energyRelative = energy/sum(energy);
    DEeach(i) =   sum(energyRelative.*myattenuation);
    meanDE  = mean(myattenuation);
    minDE = min(myattenuation);
    
    DEeach(i) = meanDE;
end

%% Make mSocketConnection with DAQ Computer

disp('Waiting for msocket communication')
%then wait for a handshake
srvsock = mslisten(42144);
masterSocket = msaccept(srvsock);
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


%% Get to the right plane
disp('Using Focus get to the right plane, click anywhere when done')
% function_Basler_Preview(Setup, 5);
function_BasPreview(Setup);

disp('Turn off Scan Image before continue')
%%

Bgdframe =  double(function_BasGetFrame(Setup,3)); %double(function_Basler_get_frames(Setup, 1 ));
Bgdframe = mean(Bgdframe,3);

[Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
[ Setup.SLM ] = Function_Start_SLM( Setup.SLM );




%%
good =0;
% while good==0
disp('Set the SLM to Divided Mode (~50) but leave it on external everything')
pwr = input('Type how much power you want to use in mW (suggest 75mW @ 10divided): ');


test=1;
Function_Feed_SLM( Setup.SLM, hololist(:,:,test));
pause(0.02);

mssend(masterSocket,[pwr/1000 DEeach(test) nTargetsList(test)]);
function_BasPreview(Setup);
% good = input('To Redo press 0, otherwise press something else');
% end

mssend(masterSocket,[0 1 1]);

invar='flush';
while ~isempty(invar)
    invar = msrecv(masterSocket,0.1);
end

%% Detrmine ROI
numTestHolos = max(containsTest);
clear dimx dimy

testFrames =[];

for i=1:numTestHolos
    holoID = find(containsTest==i,1);
    if isempty(holoID)
        disp(['Skipping Test Holo ' num2str(i)]);
    else
        Function_Feed_SLM( Setup.SLM, hololist(:,:,holoID));
        
        mssend(masterSocket,[pwr/1000 DEeach(holoID) nTargetsList(holoID)]);
        
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.1);
        end
        while ~isempty(invar)
            invar = msrecv(masterSocket,0.1);
        end
        disp('Capturing First Frame')
        
        frame =  double(function_BasGetFrame(Setup,3)); %double(function_Basler_get_frames(Setup, 1 ));
        frame = mean(frame,3);
        
        mssend(masterSocket,[0 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = msrecv(masterSocket,0.1);
        end
        while ~isempty(invar)
            invar = msrecv(masterSocket,0.1);
        end
        disp('Capturing First Frame')
        
        frame =  max(frame-Bgdframe,0);
        frame = imgaussfilt(frame,2);
        
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
        subplot(1,2,1);
        imagesc(frame);colorbar
        subplot(1,2,2);
        imagesc(smframe);colorbar
        title('Test Hologram'); 
        drawnow;    
        testFrames(:,:,i)=smframe;

    end
end

disp('Done with Test Holograms')

%%
testID=1;
Vals=zeros([1 numel(nTargetsList)]);
smallFrames = zeros([size(testFrames,1) size(testFrames,2) numel(nTargetsList)]);
randList = randperm(numel(nTargetsList));
for k=1:numel(nTargetsList)
    i=randList(k);
    
    t=tic;
    if containsTest(i)~=0
        testID = containsTest(i);
    end
    
    Function_Feed_SLM( Setup.SLM, hololist(:,:,i));
    
    mssend(masterSocket,[pwr/1000 DEeach(i) nTargetsList(i)]);
    
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.1);
    end
    % while ~isempty(invar)
    %     invar = msrecv(masterSocket,0.1);
    % end
    fprintf(['Capturing Holo ' num2str(i) '. Target ' num2str(testID)])
    
    frame =  double(function_BasGetFrame(Setup,3)); %double(function_Basler_get_frames(Setup, 1 ));
    frame = mean(frame,3);
    
    mssend(masterSocket,[0/1000 DEeach(i) nTargetsList(i)]);
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.1);
    end
    % while ~isempty(invar)
    %     invar = msrecv(masterSocket,0.1);
    % end
    
    
    frame =  max(frame-Bgdframe,0);
    frame = imgaussfilt(frame,2);
    tempx = dimx(:,testID);
    tempx(isnan(tempx))=[];
    tempy = dimy(:,testID);
    tempy(isnan(tempy))=[];
    
    smframe = frame(tempx,tempy);
    figure(1);
    subplot(1,2,1);
    imagesc(frame);colorbar
    subplot(1,2,2);
    imagesc(smframe);colorbar
    drawnow;
    
    Vals(i)=mean(smframe(:));
    smallFrames(:,:,i) = smframe; 
    
    fprintf(['. Took ' num2str(toc(t)) 's.\n']);
end

mssend(masterSocket,[0/1000 DEeach(i) nTargetsList(i)]);
disp('Done')


%%
mssend(masterSocket,'end');

%%
figure(3);clf
plot(Vals,'o')

holo = 1:numel(Vals);

scatter3(holo,Vals,nTargetsList)
ylabel('Vals')
xlabel('HoloNumber')
zlabel('NumHolos')
view(2)

out.Vals=Vals;
out.nTargetsList=nTargetsList;
out.nDistractHolos = nDistractHolos;
out.distractIter=distractIter;
out.nTestHoloToUse=nTestHoloToUse;

%%
save(['MultiTargetHolosResult_' date]);
%%

% narrowVals = squeeze(mean(mean(smallFrames(10:22,10:22,:),1),2))';
narrowVals = squeeze(mean(mean(smallFrames(13:18,13:18,:),1),2))';

%% Analysis
ValsToUse =narrowVals; Vals;
ValsToUse(containsTest==0)=[];
nTargetsListToUse =nTargetsList;
nTargetsListToUse (containsTest==0)=[];

tVals = ValsToUse(nTargetsListToUse==1);
%nDistractHolos
% t5Vals = reshape(Vals(nTargetsList==6),[distractIter,nTestHoloToUse]);
% t5Vals = reshape(Vals(nTargetsList==6),[distractIter,nTestHoloToUse]);
% t5Vals = reshape(Vals(nTargetsList==6),[distractIter,nTestHoloToUse]);
clear tDistVals
for i=1:numel(nDistractHolos)
tDistVals{i} = reshape(ValsToUse(nTargetsListToUse==nDistractHolos(i)+1),[distractIter,nTestHoloToUse]);
end



tDistMean = cell2mat(cellfun(@(x) mean(x),tDistVals,'UniformOutput',0)');
tDistSTD = cell2mat(cellfun(@(x) std(x),tDistVals,'UniformOutput',0)');
[zeros([1 size(tDistSTD,2)]); tDistSTD]
dat = [tVals; tDistMean];
figure(15);
plot([0 nDistractHolos],dat,'o-')

figure(105);clf;hold on
for i=1:size(tDistSTD,2)
    errorbar([0 nDistractHolos]+0.01*i,dat(:,i),[0; tDistSTD(:,i)], 'o-','lineWidth',2)
end

figure(115);clf;hold on
for i=1:size(tDistSTD,2)
    errorbar([0 nDistractHolos]+0.01*i,dat(:,i)./tVals(i),[0; tDistSTD(:,i)]./tVals(i), 'o-','lineWidth',2)
end
ylabel('Fluorescence (norm.)')
xlabel('Number of distractor targets')


percentResp = tDistMean./tVals;

