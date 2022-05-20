function avg = thorAcquireN(framesToAcquire, Setup)

tlCamera = Setup.cam;
% global TLCAM
% tlCamera = TLCAM;
% check to make sure cam is armed
if ~tlCamera.IsArmed
    tlCamera.Arm;
    disp('Armed camera')
end

% issue trigger
tlCamera.IssueSoftwareTrigger;
% tlCamera.ExposureTime_us=1000000;

maxPixelIntensity = double(2^tlCamera.BitDepth - 1);
numberOfFramesToAcquire = framesToAcquire;
frameCount = 0;
avgData = zeros(1080,1920,framesToAcquire); % preallocate data to average;

while frameCount < numberOfFramesToAcquire
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
            avgData(:,:,frameCount)=imageData2D';
        end
        % Release the image frame
        delete(imageFrame);
    end
end

avg = mean(avgData,3);
% figure(999)
% clf
% imagesc(avg);
% colorbar
% maxPx = max(avg, [], 'all');
% disp(['Max px was ' num2str(maxPx) ' or ' num2str(maxPx/maxPixelIntensity*100) '% of 16 bit'])

