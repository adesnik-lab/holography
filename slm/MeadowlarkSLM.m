classdef MeadowlarkSLM < handle
    properties
        is_onek = 1;
        is_loaded = 0;
        Nx = 1024;
        Ny = 1024;
        psSLM = 17e-6;
        external_pulse = 1;
        use_GPU = 1;
        max_transients = 10;
        wait_for_trigger = 0;
        ready = 0;
        timeout_ms = 5000;
    end

    properties (Constant)
        bit_depth = 12;
        is_nematic_type = 1;
        RAM_write_enable = 1;
        lut_file = 'C:\Users\Holography\Desktop\meadowlark\LUT Files\slm6257_at1035_0th_order.lut';
        lib_path = 'C:/Users/holography/Desktop/meadowlark/SDK/Blink_C_wrapper.dll';
        lib_path_h = 'C:/Users/holography/Desktop/meadowlark/SDK/Blink_C_wrapper.h';
    end

    methods

        function start(obj)

            % init C libraries
            if ~libisloaded('Blink_C_wrapper')
                fprintf('Loading SLM wrapper... ')
                loadlibrary(obj.lib_path, obj.lib_path_h);
                fprintf('done.\r')
            end

            constructed_okay = libpointer('int32Ptr', 0);
            num_boards_found = libpointer('uint32Ptr', 0);
            reg_lut = libpointer('string');

            % open SDK
            calllib('Blink_C_wrapper', 'Create_SDK', obj.bit_depth, num_boards_found, constructed_okay, obj.is_nematic_type, obj.RAM_write_enable, obj.use_GPU, obj.max_transients, reg_lut);
            
            % check creation
            if constructed_okay.value ~= 1
                disp('Blink SDK was not successfully constructed!')
                disp(calllib('Blink_C_wrapper', 'Get_last_error_message'))
                calllib('Blink_C_wrapper', 'Delete_SDK');
                disp('Deleted the SDK... Maybe just try again?')
            else
                disp('Blink SDK was successfully constructed.')
                fprintf('Found %u SLM controller(s)\n', num_boards_found.value);
                % A linear LUT must be loaded to the controller for OverDrive Plus
                fprintf('Loading LUT calibration file: %s\n', obj.lut_file)
                calllib('Blink_C_wrapper', 'Load_LUT_file', 1, obj.lut_file);
            end

            % check and report triggering status
            if obj.wait_for_trigger
                disp('Triggering active.')
            else
                disp('Triggering disabled.')
            end
            
            % all done, ready to blast
            obj.ready = 1;
            disp('SLM ready to fire!')
        end

        function blank(obj)
            mask = zeros([obj.Nx obj.Ny]);
            obj.write_mask(mask);
            disp('Blanked SLM.')
        end

        function feed(obj, holo)
            obj.write_mask(holo);
        end

        function stop(obj)
            try
                calllib('Blink_C_wrapper', 'Delete_SDK');
                disp('SLM has just been successfully turned OFF')

                if libisloaded('Blink_C_wrapper')
                    unloadlibrary('Blink_C_wrapper');
                    disp('SLM command library successfully unloaded.')
                end

            catch
                disp('Either SLM is already off, or warning for other issues...')
            end
            obj.ready = 0;
        end

    end

    methods (Access = private)
        function outcome = write_mask(obj, mask)
            calllib('Blink_C_wrapper', 'Write_image', 1, mask, obj.Nx*obj.Ny, obj.wait_for_trigger, 0, 1, 0, obj.timeout_ms);
            outcome = calllib('Blink_C_wrapper', 'ImageWriteComplete', 1, obj.timeout_ms);
        end
    end
end