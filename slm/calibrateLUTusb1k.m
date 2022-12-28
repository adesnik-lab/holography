hw = instrhwinfo('visa', 'ni');
hw.ObjectConstructorName
%%
counts = 1000;
pauseBetweenReads = 3;

v = visa('ni', 'USB0::0x1313::0x8078::P0031082::INSTR');
fopen(v);
%%

% Example usage of Blink_SDK_C.dll
% Meadowlark Optics Spatial Light Modulators
% last updated: September 10 2020

% Load the DLL
% Blink_C_wrapper.dll, Blink_SDK.dll, ImageGen.dll, FreeImage.dll and wdapi1021.dll
% should all be located in the same directory as the program referencing the
% library

addpath(genpath('c:\users\holography\desktop\meadowlark'))

if ~libisloaded('Blink_C_wrapper')
    loadlibrary('Blink_C_wrapper.dll', 'Blink_C_wrapper.h');
end

% This loads the image generation functions
if ~libisloaded('ImageGen')
    loadlibrary('ImageGen.dll', 'ImageGen.h');
end

% Basic parameters for calling Create_SDK
bit_depth = 12; %for the 1920 use 12, for the small 512x512 use 8
num_boards_found = libpointer('uint32Ptr', 0);
constructed_okay = libpointer('int32Ptr', 0);
is_nematic_type = 1;
RAM_write_enable = 1;
use_GPU = 0;
max_transients = 10;
wait_For_Trigger = 0; % This feature is user-settable; use 1 for 'on' or 0 for 'off'
external_Pulse = 0;
flip_immediate = 0; % Only supported on the 1024
timeout_ms = 5000;
RGB = 0;

%Both pulse options can be false, but only one can be true. You either generate a pulse when the new image begins loading to the SLM
%or every 1.184 ms on SLM refresh boundaries, or if both are false no output pulse is generated.
OutputPulseImageFlip = 0;
OutputPulseImageRefresh = 0; %only supported on 1920x1152, FW rev 1.8. 


% - This regional LUT file is only used with Overdrive Plus, otherwise it should always be a null string
reg_lut = libpointer('string'); %This parameter is specific to the small 512 with Overdrive, do not edit

% Call the constructor
calllib('Blink_C_wrapper', 'Create_SDK', bit_depth, num_boards_found, constructed_okay, is_nematic_type, RAM_write_enable, use_GPU, max_transients, reg_lut);

% constructed okay return of 1 is success
if constructed_okay.value ~= 1  
    disp(calllib('Blink_C_wrapper', 'Get_last_error_message'));
end
%%
if num_boards_found.value > 0 
    board_number = 1;
    disp('Blink SDK was successfully constructed');
    fprintf('Found %u SLM controller(s)\n', num_boards_found.value);
    
	% set some dimensions
	height = calllib('Blink_C_wrapper', 'Get_image_height', board_number);
    width = calllib('Blink_C_wrapper', 'Get_image_width', board_number);
	depth = calllib('Blink_C_wrapper', 'Get_image_depth', board_number); %bits per pixel
	Bytes = depth/8;
    NumDataPoints = 256;
    NumRegions = 1;
	
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
     
    %allocate arrays for our images
    Image = libpointer('uint8Ptr', zeros(width*height*Bytes,1));
	
	% ***ALWAYS*** use a blank wavefront correction when calibrating a LUT
	WFC = libpointer('uint8Ptr', zeros(width*height*Bytes,1));

    % Create an array to hold measurements from the analog input (AI) board
    AI_Intensities = zeros(NumDataPoints,2);
    
       %When calibrating you write a series of stripes to the SLM. Use the
    %summary below to set the reference, the variable grayscale, and the
    %value you step the variable grayscale by.
    
    %1920x1152, reference = 0, variable = 0 to 255 in steps of +1
    %512x512 8-bit, reference = 255, variable = 255 to 0 in steps of -1
    %512x512 16-bit, reference = 63353, variable = 63353 to 0 in steps of -256
    %1024x1024, reference = 255, variable = 255 to 0 in steps of -1
    Reference = 255;
    Variable = 255;
    StepBy = -1;
    
    % Use a high frequency grating to separate the 0th and 1st orders, a
    % period of 8 is generally good
    PixelsPerStripe = 8;
    testGreyLevel = 10;
    PixelValue = 0;
    Region = 0;
    calllib('ImageGen', 'Generate_Stripe', Image, WFC, width, height, depth, Reference, testGreyLevel, PixelsPerStripe, RGB);
    calllib('ImageGen', 'Mask_Image', Image, width, height, depth, Region, NumRegions, RGB);
    calllib('Blink_C_wrapper', 'Write_image', board_number, Image, width*height*Bytes, wait_For_Trigger, flip_immediate, OutputPulseImageFlip, OutputPulseImageRefresh, timeout_ms);
	
    input('Press enter when you align everything. There should be stripes on the SLM.')
    
	% begin with the SLM blank
    calllib('Blink_C_wrapper', 'Write_image', board_number, zeros([1024, 1024], 'uint8'), width*height*Bytes, wait_For_Trigger, flip_immediate, OutputPulseImageFlip, OutputPulseImageRefresh, timeout_ms);
	calllib('Blink_C_wrapper', 'ImageWriteComplete', board_number, timeout_ms);
	
    pause(1)
    
 
    %loop through each region
    for Region = 0:(NumRegions-1)
        
        figure(1)
        clf
      
        AI_Index = 1;
        %loop through each graylevel
        
        c=0;
        cs=[];
        vals = [];
        
        for TestPoint = 0:(NumDataPoints-1)
            %Generate the stripe pattern and mask out current region
            calllib('ImageGen', 'Generate_Stripe', Image, WFC, width, height, depth, Reference, Variable, PixelsPerStripe, RGB);
            calllib('ImageGen', 'Mask_Image', Image, width, height, depth, Region, NumRegions, RGB);
            
            %Step the variable grayscale
            Variable = Variable + StepBy;
            
            %write the image
            calllib('Blink_C_wrapper', 'Write_image', board_number, Image, width*height*Bytes, wait_For_Trigger, flip_immediate, OutputPulseImageFlip, OutputPulseImageRefresh, timeout_ms);
            
            %let the SLM settle for 10 ms
            pause(0.1);
            c = c+1;
            if c == 1
                input('let power meter stabilize, then press enter to start')
            end
            %YOU FILL IN HERE...FIRST: read from your specific AI board, note it might help to clean up noise to average several readings
            %SECOND: store the measurement in your AI_Intensities array
            pause(pauseBetweenReads)
            fprintf(v,['sense:average:count ',num2str(counts)]);
            set(v,'timeout',3+1.1*counts*3/1000)
            ret = query(v,'read?');
            val = str2num(ret)*1000;
            disp(['Got Measurement (' num2str(Variable) '): ' num2str(val)])

            cs = [cs c];
            vals = [vals val];
            scatter(cs, vals)
            xlim([0 255])
            
            AI_Intensities(AI_Index, 1) = TestPoint; %This is the varable graylevel you wrote to collect this data point
            AI_Intensities(AI_Index, 2) = val; % HERE YOU NEED TO REPLACE 0 with YOUR MEASURED VALUE FROM YOUR ANALOG INPUT BOARD
 
            AI_Index = AI_Index + 1;
        
        end
        
        disp('all done. saving...')
        
        usrhome = getenv('USERPROFILE');
        
        % dump the AI measurements to a csv file
        filename = ['Raw' num2str(Region) '.csv'];
        filepath = fullfile(usrhome, 'Desktop', 'meadowlark', 'LUT files Raw', filename)
        csvwrite(filepath, AI_Intensities);  
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