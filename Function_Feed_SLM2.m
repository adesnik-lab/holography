function [ outcome ] = Function_Feed_SLM2( SLM,frame)
calllib('Blink_C_wrapper', 'Write_image', 1, frame, 1920*1152, SLM.wait_For_Trigger, SLM.external_Pulse, SLM.timeout_ms);
outcome = calllib('Blink_C_wrapper', 'ImageWriteComplete', 1, SLM.timeout_ms);
end

