
function [Setup] = function_startBasCam(varargin)

if nargin<1
    Setup=[];
else
    Setup=varargin{1};
end

if isfield(Setup,'useThorCam') && Setup.useThorCam==1
    
    Setup.cam = thorSetupCam;
    Setup.cam.ExposureTime_us = 10000;
    
    if isfield(Setup,'camExposureTime')
        Setup.cam.ExposureTime_us  =Setup.camExposureTime;
    end
    if isfield(Setup,'maxFramesPerAcquire') && Setup.maxFramesPerAcquire>0
        disp(['Changing the number of frames per trigger to: ' num2str(Setup.maxFramesPerAcquire)]);
        Setup.cam.Disarm;
        Setup.cam.MaximumNumberOfFramesToQueue = Setup.maxFramesPerAcquire;
        Setup.cam.FramesPerTrigger_zeroForUnlimited = Setup.maxFramesPerAcquire;
        Setup.cam.Arm;
    end
else
    
    if isfield(Setup,'cam')
        vid = Setup.cam;
    else
        %     vid = videoinput('winvideo', 1, 'Y800_1024x768');
        %     vid = videoinput('winvideo', 1, 'Y800_2048x1536');
        
        %vid = videoinput('winvideo', 1, 'Y800_2592x1944');
        try
            vid = videoinput('winvideo', 1, 'Y800_1280x1024');
            % vid = videoinput('winvideo', 1, 'Y800_1024x768');
        catch
            vid = videoinput('winvideo', 1, 'Y800_2048x1536');
        end
        
        src = getselectedsource(vid);
        vid.ReturnedColorspace = 'grayscale';
        vid.FramesPerTrigger = 1;
        triggerconfig(vid, 'manual');
        vid.TriggerRepeat = Inf;
        src.Exposure = -3; %sensitivity of camera range -14 to 0. older camera is 0. -Ian 3/6/20
        src.Gain = 500; % wasn't set previously, value was 136, 7/19/20 WH
    end
    
    start(vid);
    
    Setup.cam = vid;
end