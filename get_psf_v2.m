%% Pathing and setup

clear
close all
clc

tBegin = tic;

disp('Setting up stuff...');

makePaths()

Setup = function_loadparameters3();
Setup.CGHMethod = 2;
Setup.GSoffset = 0;
Setup.verbose = 0;
Setup.useGPU = 1;

if Setup.useGPU
    fprintf('Getting gpu... '); %this can sometimes take a while at initialization
    g = gpuDevice;
    fprintf('done.\r')
end

slm = MeadowlarkSLM();
slm.start()

sutter = sutterController();

bas = bascam();
bas.start()

disp('Setup complete.')
%% look for objective in 1p

bas.preview()

%% connect to DAQ computer

%run this first then code on daq

laser = LaserSocket(42133);


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

pwr = 2.4;
slmCoords = [.55 .45 0 1];


disp(['Individual hologram power set to ' num2str(pwr) 'mW.'])

DEestimate = DEfromSLMCoords(slmCoords);
disp(['Diffraction Estimate for this spot is: ' num2str(DEestimate)])

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
slm.feed(Holo)
laser.set_power(pwr)
bas.preview()
laser.set_power(0)

%% check pixel val
bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

laser.set_power(pwr)
data = bas.grab(10);
laser.set_power(0)

data = mean(data,3);
frame = max(data-bgd, 0);
[x,y] = function_findcenter(frame);

frame_crop = data(x-20:x+20, y-20:y+20);

figure(6)
clf
imagesc(frame_crop)
axis square
colorbar
% clim([0 255])
xL=xlim;
yL=ylim;
mx = max(frame_crop,[],"all");
str = ['Peak Intensity: ' num2str(mx) ' A.U.'];
text(0.03*xL(2),0.03*yL(2),str,'HorizontalAlignment','left','VerticalAlignment','top', 'Color','w')

%% Measure background noise in camera


bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

meanBgd = mean(bgd_frames, 'all');
stdBgd = std(double(bgd_frames), [], 'all');

%%Collect PSF
nframesCapture = 5;

figure(1); clf

sutter.setRef()

sz = size(bgd);
UZ = linspace(-50,50,31);

dataUZ = zeros([sz numel(UZ)]);

disp('Collecting PSF.')

for i=1:numel(UZ)

    sutter.moveZ(UZ(i))

    if i==1
        pause(1)
    else
        pause(0.1);
    end

    laser.set_power(pwr)
    data = bas.grab(nframesCapture);
    laser.set_power(0)

     % process frame grabs
    data = mean(data, 3);
    frame = max(data-bgd, 0);
    frame = imgaussfilt(frame, 2);
    dataUZ(:,:,i) =  frame;

    % plot live data
    figure(1)
    subplot(1,2,1)
    imagesc(frame);
    title(['Frame ' num2str(i)])
    colorbar
%     axis equal
%     clim([0 255])

    subplot(1,2,2)
    imagesc(max(dataUZ,[],3));
    colorbar
%     clim([0 255])
    title('Max Projection')
%     axis equal
    drawnow
end

sutter.moveToRef()

disp('Done collecting PSF stack.')

%%Calc FWHM

mxProj = max(dataUZ,[],3);
[ x,y ] =function_findcenter(mxProj);

range = 2;

dimx = max((x-range),1):min((x+range),size(frame,1));
dimy = max((y-range),1):min((y+range),size(frame,2));

thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
[a, peakPlane ] = max(thisStack);
peakFrame = dataUZ(:,:,peakPlane);
thisStack = thisStack - min(thisStack);

% determine pxPerMu with sutter
muUsed = 50;

disp('Determining pxPerMu.')

xyz = [0 0 UZ(peakPlane)];
sutter.moveTo(xyz)
pause(1)

laser.set_power(pwr)
data = bas.grab(nframesCapture);
data = mean(data, 3);
laser.set_power(0)

frame = max(data-bgd, 0);
pos1 = imgaussfilt(frame, 2);

xyz = [0 muUsed UZ(peakPlane)];
sutter.moveTo(xyz)
pause(1)

laser.set_power(pwr)
data = bas.grab(nframesCapture);
data = mean(data, 3);
laser.set_power(0)

frame = max(data-bgd, 0);
pos2 = imgaussfilt(frame, 2);

[ x1,y1 ] =function_findcenter(pos1 );
[ x2,y2 ] =function_findcenter(pos2 );
distance = pdist([x1 y1; x2 y2]);
pxPerMu = distance/muUsed;
disp(['Determined to be ' num2str(pxPerMu) ' px per mu'])

% figure(4)
% clf
% subplot(1,2,1)
% imagesc(pos1)
% subplot(1,2,2)
% imagesc(pos2)

sutter.moveToRef()

xLine = double(peakFrame(x,:));
yLine = double(peakFrame(:,y))';
[row,col] = size(frame);
xSize = linspace(1,col,numel(xLine))./pxPerMu; 
ySize = linspace(1,row,numel(yLine))./pxPerMu; 

f1 = fit(xSize', xLine', 'gauss1');
xValue =f1.a1;
xDepth =f1.b1;
xFWHM = 2*sqrt(2*log(2))*f1.c1/sqrt(2);

ff = fit(UZ', thisStack, 'gauss1');
peakValue =ff.a1;
peakDepth =ff.b1;
peakFWHM = 2*sqrt(2*log(2))*ff.c1/sqrt(2);

mimg = mean(dataUZ,3);

% caclulations
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

disp(['Max camera val: ' num2str(max(peakFrame(:)))])


