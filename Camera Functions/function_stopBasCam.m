function [Setup] = function_stopBasCam(Setup)

if isfield(Setup,'useThorCam') && Setup.useThorCam==1
    thorStop(Setup.cam)
    Setup = rmfield(Setup,'cam');
else
    vid = Setup.cam;
    
    stop(vid);
    delete(vid);
    Setup = rmfield(Setup,'cam');
end