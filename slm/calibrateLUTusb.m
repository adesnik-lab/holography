hw = instrhwinfo('visa', 'ni');
hw.ObjectConstructorName
%%
counts = 500;

v = visa('ni', 'USB0::0x1313::0x8078::P0031082::INSTR');
fopen(v);
%%
% Run this first, then run the DiffractiveLUT.exe
% This is the manual version, would be easy to make automatic with USB
% cable for power meter.

directions = ['\n\nTo collect the LUT using the zero order, place a pinhole at \n'...
              'the zero order (remove zero-order) and close it so it is small \n'...
              'enough to only let the light from the zero-order in. Then place the power \n'...
              'meter behind the pinhole.\n'...
              '\nThen collect the powers for each. Units (mW or W) don''t matter, just \n'...
              'be consistent. Use enough power to get a good curve but not enough to burn the SLM.\n'];

fprintf(directions)
input('Press ENTER to continue: ')

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
    

    
    %set some dimensions
    height = calllib('Blink_C_wrapper', 'Get_image_height', board_number);
    width = calllib('Blink_C_wrapper', 'Get_image_width', board_number);
    NumDataPoints = 256;
    NumRegions = 1;
    
    % To measure the raw response we want to disable the LUT by loading a linear LUT
    % To measure the raw optical response we want to linearly increment the voltage on the pixels by using a linear LUT
    if ((width == 512) && (depth == 8))
		calllib('Blink_C_wrapper', 'Load_LUT_file', board_number, 'C:\\Program Files\\Meadowlark Optics\\Blink OverDrive Plus\\LUT Files\\512x512_linearVoltage.LUT');
    end
    if ((width == 512) && (depth == 16))
		calllib('Blink_C_wrapper', 'Load_LUT_file', board_number, 'C:\\Program Files\\Meadowlark Optics\\Blink OverDrive Plus\\LUT Files\\512x512_16bit_linearVoltage.LUT');
    end
    if width == 1920
		calllib('Blink_C_wrapper', 'Load_LUT_file', board_number, 'C:\\Program Files\\Meadowlark Optics\\Blink OverDrive Plus\\LUT Files\\1920x1152_linearVoltage.LUT');
    end
    if width == 1024
		calllib('Blink_C_wrapper', 'Load_LUT_file', board_number, 'C:\\Program Files\\Meadowlark Optics\\Blink OverDrive Plus\\LUT Files\\1024x1024_linearVoltage.LUT');
    end
    calllib('Blink_C_wrapper', 'Load_LUT_file',board_number, lut_file);
    
    %allocate arrays for our images
    Image = libpointer('uint8Ptr', zeros(width*height,1));

    % Create an array to hold measurements from the analog input (AI) board
    AI_Intensities = zeros(NumDataPoints,2);
    
    % ***ALWAYS*** use a blank wavefront correction when calibrating a LUT
	WFC = libpointer('uint8Ptr', zeros(width*height*Bytes,1));
    
    % Generate a blank wavefront correction image, you should load your
    % custom wavefront correction that was shipped with your SLM.
    PixelsPerStripe = 8;
    PixelValue = 0;
    Region=0;
    calllib('ImageGen', 'Generate_Solid', Image, width, height, PixelValue);
    calllib('Blink_C_wrapper', 'Write_image', board_number, Image, width*height, wait_For_Trigger, external_Pulse, timeout_ms);
	calllib('Blink_C_wrapper', 'ImageWriteComplete', board_number, timeout_ms);
    
    testGreyLevel = 10;
    calllib('ImageGen', 'Generate_Stripe', Image, width, height, PixelValue, testGreyLevel, PixelsPerStripe);
    calllib('ImageGen', 'Mask_Image', Image, width, height, Region, NumRegions);
    calllib('Blink_C_wrapper', 'Write_image', board_number, Image, width*height, wait_For_Trigger, external_Pulse, timeout_ms);
	
    %%
    
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
            
            %let the SLM settle for 100 ms
            pause(0.1);
            
            c = c+1;
            %YOU FILL IN HERE...FIRST: read from your specific AI board, note it might help to clean up noise to average several readings
            %SECOND: store the measurement in your AI_Intensities array
            AI_Intensities(AI_Index, 1) = Gray; %This is the varable graylevel you wrote to collect this data point
            
            pause(3)
            fprintf(v,['sense:average:count ',num2str(counts)]);
            set(v,'timeout',3+1.1*counts*3/1000)
            ret = query(v,'read?');
            val = str2num(ret)*1000;
            disp(['Got Measurement (' num2str(Gray) '): ' num2str(val)])

            AI_Intensities(AI_Index, 2) = val;
            cs = [cs c];
            vals = [vals val];
            scatter(cs, vals)
            xlim([0 255])

            AI_Index = AI_Index + 1;
        
        end
        disp('all done. saving...')
        
        % dump the AI measurements to a csv file and mat file
        % you only need the csv but it's helpful to have the mat in case
        % you need to edit it
        
        if ispc
            usrhome = getenv('USERPROFILE');
        else
            usrhome = getenv('HOME');
        end
        
        lut_folder = fullfile(usrhome, 'Desktop', 'SLM_LUT_files');
        status = mkdir(lut_folder);
        
        
        filename_csv = fullfile(lut_folder, ['Raw' num2str(Region) '.csv']);
        csvwrite(filename_csv, AI_Intensities);
        
        filename_mat = fullfile(lut_folder, 'RawValuesLUT.mat');
        save(filename_mat,'AI_Intensities')
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

disp('data saved!')