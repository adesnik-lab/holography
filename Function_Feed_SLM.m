function [ outcome ] = Function_Feed_SLM( SLM,frame)
if SLM.is_onek
    calllib('Blink_C_wrapper', 'Write_image', 1, frame, SLM.Nx*SLM.Ny, SLM.wait_For_Trigger,0, 1, 0, SLM.timeout_ms);
else
    calllib('Blink_C_wrapper', 'Write_image', 1, frame, 1920*1152, SLM.wait_For_Trigger, SLM.external_Pulse, SLM.timeout_ms);
end
outcome = calllib('Blink_C_wrapper', 'ImageWriteComplete', 1, SLM.timeout_ms);
end

