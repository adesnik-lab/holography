function [frames, t] = function_BasGetFrame(Setup,N)
if isfield(Setup,'useThorCam') && Setup.useThorCam==1
 tic;
 if Setup.cam.FramesPerTrigger_zeroForUnlimited ~= N
     disp(['Changing the number of frames per trigger to: ' num2str(N)]);
     Setup.cam.Disarm;
     Setup.cam.MaximumNumberOfFramesToQueue = N;
     Setup.cam.FramesPerTrigger_zeroForUnlimited = N;
     Setup.cam.Arm;
 end
 
 timeToPause = N * double(Setup.cam.ExposureTime_us) / 1e6;
 pause(timeToPause);
 
 frames = thorAcquireN(N, Setup);
 t = toc;
else

N=N+1; %throwout first frame

if isfield(Setup,'cam')
    vid = Setup.cam;
else
    disp('Problem with setup')
    return
% %     vid = videoinput('winvideo', 1, 'Y800_1024x768');
%     vid = videoinput('winvideo', 1, 'Y800_2048x1536');
%     src = getselectedsource(vid);
%     vid.ReturnedColorspace = 'grayscale';
end

%clear buffer

while vid.FramesAvailable~=0
    flushdata(vid);
end

stopFlag=0;
if strcmp(vid.Running,'on')
    if vid.FramesPerTrigger~=N ||...
            vid.TriggerRepeat ~= Inf ||...
            ~strcmp(vid.TriggerType,'manual')
        
        stopFlag =1;
    end
else
    stopFlag = 1;
end
if stopFlag
    disp('Changing Parameters...')
    stop(vid);
    vid.FramesPerTrigger = N;
    triggerconfig(vid, 'manual');
    vid.TriggerRepeat = Inf;
    
    start(vid);
end
tic
trigger(vid);
x = getdata(vid);
t = toc;
% size(x)

frames = squeeze(x(:,:,1,2:end));
end