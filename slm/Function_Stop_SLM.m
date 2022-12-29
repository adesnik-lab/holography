function [SLM ] = Function_Stop_SLM( SLM )
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

SLM.State = 0;

end

