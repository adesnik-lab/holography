function thorStop(tlCamera,other)
    global TLCAMSDK TLCAM SERIALNUMBERS

if nargin<1
    tlCamera = TLCAM;
    tlCameraSDK = TLCAMSDK;
    serialNumbers = SERIALNUMBERS;
elseif nargin<2
     serialNumbers = SERIALNUMBERS;
     tlCameraSDK = TLCAMSDK;
else
    serialNumbers = other.serialNumbers;
    tlCameraSDK = other.tlCameraSDK;
end

% Stop continuous image acquisition
disp('Stopping continuous image acquisition.');
tlCamera.Disarm;

% Release the TLCamera
disp('Releasing the camera');
tlCamera.Dispose;
delete(tlCamera);


delete(serialNumbers);

% Release the TLCameraSDK.
tlCameraSDK.Dispose;
delete(tlCameraSDK);