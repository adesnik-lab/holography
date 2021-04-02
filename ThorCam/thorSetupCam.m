function tlCamera = thorSetupCam

% Adapted from Thorlabs
try
    thorStop;
    disp('stopped running cam')
catch
end
% clear
% close all

% Load TLCamera DotNet assembly. The assembly .dll is assumed to be in the
% same folder as the scripts.
oldLoc = cd; 
cam_path = 'C:\Users\Holography\Desktop\holography\ThorCam\';
pause(0.1)
cd(cam_path);
NET.addAssembly([cam_path, 'Thorlabs.TSI.TLCamera.dll']);
% NET.addAssembly('C:\Users\ian\Dropbox\code\Scientific Camera Interfaces\Thorlabs.TSI.TLCamera.dll');

disp('Dot NET assembly loaded.');

tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;

% Get serial numbers of connected TLCameras.
serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
disp([num2str(serialNumbers.Count), ' camera was discovered.']);

disp('Opening the first camera')
tlCamera = tlCameraSDK.OpenCamera(serialNumbers.Item(0), false);

% Set exposure time and gain of the camera.
tlCamera.ExposureTime_us = 1000000;

% Check if the camera supports setting "Gain"
gainRange = tlCamera.GainRange;
if (gainRange.Maximum > 0)
    tlCamera.Gain = 0;
end

% Set the FIFO frame buffer size. Default size is 1.
tlCamera.MaximumNumberOfFramesToQueue = 5;
% arm the camera
tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
tlCamera.Arm;
disp('Camera armed.')

other.serialNumbers = serialNumbers;
other.tlCameraSDK = tlCameraSDK;

global TLCAMSDK TLCAM SERIALNUMBERS
TLCAMSDK = tlCameraSDK;
TLCAM = tlCamera;
SERIALNUMBERS= serialNumbers;
cd(oldLoc)


