%% exclude trials

slmXYZ = slmXYZBackup2;
basXYZ = basXYZBackup2;
basVal = basValBackup2;
excludeTrials = all(basXYZ(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

excludeTrials = excludeTrials | basVal>260; %max of this camera is 255

basDimensions = size(Bgdframe);
excludeTrials = excludeTrials | basXYZ(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | any(isnan(basXYZ(:,:)));
%excludeTrials = excludeTrials | basVal<5; %8/3 hayley add

slmXYZBackup = slmXYZ(:,~excludeTrials);
basXYZBackup = basXYZ(:,~excludeTrials);
basValBackup = basVal(:,~excludeTrials);
%% fit SLM to Camera
%use model terms

basXYZ = basXYZBackup;
slmXYZ = slmXYZBackup;
basVal = basValBackup;

disp('Fitting SLM to Camera')
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1 ; 1 1 1 ;...
    2 0 0; 0 2 0; 0 0 2;  ...
    2 0 1; 2 1 0; 0 2 1; 1 2 0; 0 1 2;  1 0 2; ... ];  %XY spatial calibration model for Power interpolations
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2;];
reOrder = randperm(size(slmXYZ,2));
slmXYZ = slmXYZ(:,reOrder);
basXYZ = basXYZ(:,reOrder);

holdback = 100;%50;

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = (basXYZ(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

errScalar = 2.5;%2.5;
[SLMtoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);

Ask = refAsk;
True = refGet;
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(103);clf
subplot(1,2,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
% legend('Measured targets', 'Estimated Targets');
title({'Reference Data'; 'SLM to Camera'})

refRMS = sqrt(sum((Get-True).^2,2));
subplot(1,2,2)
scatter3(True(:,1),True(:,2),True(:,3),[],refRMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title({'Reference Data'; 'RMS Error in position'})
caxis([0 30])


Ask = (slmXYZ(1:3,end-holdback:end))';
True = (basXYZ(1:3,end-holdback:end))';
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(101);clf
subplot(1,3,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')


ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
legend('Measured targets', 'Estimated Targets');
title('SLM to Camera')

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp('Error based on Holdback Data...')
disp(['The RMS error: ' num2str(meanRMS) ' pixels for SLM to Camera']);

pxPerMu = size(frame,1) / 1000; %really rough approximate of imaging size

disp(['Thats approx ' num2str(meanRMS/pxPerMu) ' um']);

xErr = sqrt(sum((Get(:,1)-True(:,1)).^2,2));
yErr = sqrt(sum((Get(:,2)-True(:,2)).^2,2));
zErr = sqrt(sum((Get(:,3)-True(:,3)).^2,2));

disp('Mean:')
disp(['X: ' num2str(mean(xErr)/pxPerMu) 'um. Y: ' num2str(mean(yErr)/pxPerMu) 'um. Z: ' num2str(mean(zErr)) 'um.']);
disp('Max:')
disp(['X: ' num2str(max(xErr)/pxPerMu) 'um. Y: ' num2str(max(yErr)/pxPerMu) 'um. Z: ' num2str(max(zErr)) 'um.']);

subplot(1,3,2)
scatter3(True(:,1),True(:,2),True(:,3),[],RMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title('RMS Error in position')


refAsk = (basXYZ(1:3,1:end-holdback))';
refGet = (slmXYZ(1:3,1:end-holdback))';

%  camToSLM = function_3DCoC(refAsk,refGet,modelterms);

[camToSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);


Ask = (basXYZ(1:3,end-holdback:end))';
True = (slmXYZ(1:3,end-holdback:end))';
Get = function_Eval3DCoC(camToSLM,Ask);


figure(101);
subplot(1,3,3)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Depth units')
legend('Measured targets', 'Estimated Targets');
title('Camera to SLM')

% RMS = sqrt(sum((Get-True).^2,2));
% meanRMS = nanmean(RMS);
% 
% disp(['The RMS error: ' num2str(meanRMS) ' SLM units for Camera to SLM']);




CoC.camToSLM=camToSLM;
CoC.SLMtoCam = SLMtoCam;

out.CoC=CoC;
out.CoCmodelterms = modelterms;

rtXYZ = function_Eval3DCoC(SLMtoCam,function_Eval3DCoC(camToSLM,basXYZ(1:3,end-holdback:end)'));

err = sqrt(sum((rtXYZ - basXYZ(1:3,end-holdback:end)').^2,2));
meanRTerr = nanmean(err);
disp(['The Mean Round Trip RMS error: ' num2str(meanRTerr) ' pixels (' num2str(meanRTerr/pxPerMu) ' um) camera to SLM to camera']);


disp('Fitting Power as a function of SLM')
%  modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%      1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%      2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;   ];  %XY spatial calibration model for Power interpolations
slmXYZ = slmXYZBackup;
basVal = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

intVal = basVal;
intVal = sqrt(intVal); %convert fluorescence intensity (2P) to 1P illumination intensity
intVal=intVal./max(intVal(:));

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = intVal(1:end-holdback);

SLMtoPower =  polyfitn(refAsk,refGet,modelterms);

Ask = (slmXYZ(1:3,end-holdback:end))';
True = intVal(end-holdback:end)';

Get = polyvaln(SLMtoPower,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);

figure(1);clf
subplot(2,3,1)
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],intVal,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Measured Power (converted to 1p)')
colorbar
axis square

subplot(2,3,2)
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)'),'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated Power Norm.')
colorbar
axis square

subplot(2,3,4)
% scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal','filled');
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],basVal,'filled');

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Raw Fluorescence')
colorbar
axis square

subplot(2,3,5)
c = sqrt((polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal').^2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error RMS (A.U.)')
colorbar
axis square

subplot(2,3,3)
c = (polyvaln(SLMtoPower,slmXYZ(1:3,:)').^2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated 2P Power')
colorbar
axis square

subplot(2,3,6)
normVal = basVal./max(basVal(:));

c = (polyvaln(SLMtoPower,slmXYZ(1:3,:)').^2)-normVal';
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error 2P Power')
colorbar
axis square


disp(['The RMS error: ' num2str(meanRMS) ' A.U. Power Estimate']);
disp(['The Max power error: ' num2str(max(RMS)*100) '% of request']);

CoC.SLMtoPower = SLMtoPower;
out.CoC = CoC;
out.powerFitmodelTerms = modelterms;