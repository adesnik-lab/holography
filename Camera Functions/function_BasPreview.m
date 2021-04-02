function [Setup] = function_BasPreview(varargin)

if nargin<1
    Setup=[];
else
    Setup=varargin{1};
end

if isfield(Setup,'useThorCam') && Setup.useThorCam==1
    if ~isfield(Setup,'cam')
        [Setup.cam, other] = thorSetupCam;        
        thorPreview(Setup)
        thorStop(Setup.cam)
        
    else
        
        if Setup.cam.FramesPerTrigger_zeroForUnlimited~=0
            OldFramesPerTrigger = Setup.cam.FramesPerTrigger_zeroForUnlimited;
            Setup.cam.Disarm;
            Setup.cam.FramesPerTrigger_zeroForUnlimited = 0;
            Setup.cam.Arm;
            
            thorPreview(Setup)
            
            Setup.cam.Disarm;
            Setup.cam.FramesPerTrigger_zeroForUnlimited = OldFramesPerTrigger;
            Setup.cam.Arm;
        else
            thorPreview(Setup)
        end
        
        Setup.cam.ExposureTime_us = Setup.camExposureTime;
        
    end
else
    if isfield(Setup,'cam')
        vid = Setup.cam;
    else
        %     vid = videoinput('winvideo', 1, 'Y800_2048x1536');
        vid = videoinput('winvideo', 1, 'Y800_1024x768');
        %     vid = videoinput('winvideo', 1, 'Y800_800x600');
        
        src = getselectedsource(vid);
        vid.ReturnedColorspace = 'grayscale';
    end
    
    Setup.cam = vid;
    h = preview(vid);
    %
    % while isgraphics(h)
    %     continue
    % end
    %     stoppreview(vid);
    %     closepreview;
    %
    dummy = input('Press Enter to Continue');
    
    
    stoppreview(vid);
    closepreview;
end