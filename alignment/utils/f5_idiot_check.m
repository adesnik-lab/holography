function f5_idiot_check()
% guards against hitting F5 instead of F9

in = input('Run calibration? ..or did you hit F5 by mistake? (y/N) ','s');
if ~strcmp(in,'y')
    error('ERROR- Not starting calibration.')
else
    disp('Beginning calibration... have fun!')
end