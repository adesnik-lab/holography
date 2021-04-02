function thorPreview(Setup)
% if nargin<1
%     [tlCamera, other] = thorSetupCam;
%     unMountWhenDone = 1;
% else
%     unMountWhenDone = 0;
% end

tlCamera = Setup.cam;

figHandle = figure('units','normalized','outerposition',[0 0 1 1]);

% Start continuous image acquisition
disp('Starting continuous image acquisition.');
%tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
%tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
% try
%     tlCamera.Disarm;
% end
% tlCamera.Arm;
tlCamera.IssueSoftwareTrigger;
maxPixelIntensity = double(2^tlCamera.BitDepth - 1);
tlCamera.ExposureTime_us=200000;%25000;
numberOfFramesToAcquire = 1;
frameCount = 0;

while isgraphics(figHandle)
    % Check if image buffer has been filled
    if (tlCamera.NumberOfQueuedFrames > 0)
        
        % Get the pending image frame.
        imageFrame = tlCamera.GetPendingFrameOrNull;
        
        if ~isempty(imageFrame)
            frameCount = frameCount + 1;
            
            % Get the image data as 1D uint16 array
            imageData = uint16(imageFrame.ImageData.ImageData_monoOrBGR);            
            imageHeight = imageFrame.ImageData.Height_pixels;
            imageWidth = imageFrame.ImageData.Width_pixels;
            imageData2D = reshape(imageData, [imageWidth, imageHeight]);
%             figure(figHandle)
            imagesc(imageData2D')
            colormap(gray)
            colorbar
        end
        
        % Release the image frame
        delete(imageFrame);
    end
    drawnow;
    
end

% % Stop continuous image acquisition
disp('Stopping continuous image acquisition.');
% tlCamera.Disarm;
% 
% if unMountWhenDone
%     stopCamera(tlCamera,other)
% end