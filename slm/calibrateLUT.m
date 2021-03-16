% Run this first, then run the DiffractiveLUT.exe

% Example usage of Blink_SDK_C.dll
% Meadowlark Optics Spatial Light Modulators
% last updated: April 6 2018

% Load the DLL
% Blink_C_wrapper.dll, Blink_SDK.dll, ImageGen.dll, FreeImage.dll and wdapi1021.dll
% should all be located in the same directory as the program referencing the
% library
if ~libisloaded('Blink_C_wrapper')
    loadlibrary('Blink_C_wrapper.dll', 'Blink_C_wrapper.h');
end

% This loads the image generation functions
if ~libisloaded('ImageGen')
    loadlibrary('ImageGen.dll', 'ImageGen.h');
end

% Basic parameters for calling Create_SDK
bit_depth = 12;
num_boards_found = libpointer('uint32Ptr', 0);
constructed_okay = libpointer('int32Ptr', 0);
is_nematic_type = 1;
RAM_write_enable = 1;
use_GPU = 0;
max_transients = 10;
wait_For_Trigger = 0; % This feature is user-settable; use 1 for 'on' or 0 for 'off'
external_Pulse = 0;
timeout_ms = 5000;
reg_lut = libpointer('string');

% Call the constructor
calllib('Blink_C_wrapper', 'Create_SDK', bit_depth, num_boards_found, constructed_okay, is_nematic_type, RAM_write_enable, use_GPU, max_transients, reg_lut);

% Convention follows that of C function return values: 0 is success, nonzero integer is an error
if constructed_okay.value ~= 0  
    disp('Blink SDK was not successfully constructed');
    disp(calllib('Blink_C_wrapper', 'Get_last_error_message'));
    calllib('Blink_C_wrapper', 'Delete_SDK');
else
    board_number = 1;
    disp('Blink SDK was successfully constructed');
    fprintf('Found %u SLM controller(s)\n', num_boards_found.value);
    
    % To measure the raw response we want to disable the LUT by loading a linear LUT
    lut_file = 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\linear.LUT';
    calllib('Blink_C_wrapper', 'Load_LUT_file',board_number, lut_file);
    
    %set some dimensions
    height = calllib('Blink_C_wrapper', 'Get_image_height', board_number);
    width = calllib('Blink_C_wrapper', 'Get_image_width', board_number);
    NumDataPoints = 256;
    NumRegions = 1;
    
    %allocate arrays for our images
    Image = libpointer('uint8Ptr', zeros(width*height,1));

    % Create an array to hold measurements from the analog input (AI) board
    AI_Intensities = zeros(NumDataPoints,2);
    
    % Generate a blank wavefront correction image, you should load your
    % custom wavefront correction that was shipped with your SLM.
    PixelValue = 0;
    calllib('ImageGen', 'Generate_Solid', Image, width, height, PixelValue);
    calllib('Blink_C_wrapper', 'Write_image', board_number, Image, width*height, wait_For_Trigger, external_Pulse, timeout_ms);
	calllib('Blink_C_wrapper', 'ImageWriteComplete', board_number, timeout_ms);
	
    PixelsPerStripe = 8;
    %loop through each region
    for Region = 0:(NumRegions-1)
      
        AI_Index = 1;
        %loop through each graylevel
        c=0;
        cs=[];
        vals = [];
        for Gray = 0:(NumDataPoints-1)
            %Generate the stripe pattern and mask out current region
            calllib('ImageGen', 'Generate_Stripe', Image, width, height, PixelValue, Gray, PixelsPerStripe);
            calllib('ImageGen', 'Mask_Image', Image, width, height, Region, NumRegions);
            
            %write the image
            calllib('Blink_C_wrapper', 'Write_image', board_number, Image, width*height, wait_For_Trigger, external_Pulse, timeout_ms);
            
            %let the SLM settle for 10 ms
            pause(0.01);
            
            c = c+1;
            %YOU FILL IN HERE...FIRST: read from your specific AI board, note it might help to clean up noise to average several readings
            %SECOND: store the measurement in your AI_Intensities array
            AI_Intensities(AI_Index, 1) = Gray; %This is the varable graylevel you wrote to collect this data point
            val = input('Enter power in mW: '); % HERE YOU NEED TO REPLACE 0 with YOUR MEASURED VALUE FROM YOUR ANALOG INPUT BOARD
            AI_Intensities(AI_Index, 2) = val;
            cs = [cs c];
            vals = [vals val];
            scatter(cs, vals)
            
            AI_Index = AI_Index + 1;
        
        end
        
        % dump the AI measurements to a csv file
        filename = ['Raw' num2str(Region) '.csv'];
        csvwrite(filename, AI_Intensities);
        save('RawValuesLUT','AI_Intensities')
    end
	
     
    % Always call Delete_SDK before exiting
    calllib('Blink_C_wrapper', 'Delete_SDK');
end

%destruct
if libisloaded('Blink_C_wrapper')
    unloadlibrary('Blink_C_wrapper');
end

if libisloaded('ImageGen')
    unloadlibrary('ImageGen');
end