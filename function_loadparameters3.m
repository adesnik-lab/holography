function Setup = function_loadparameters3()

Setup.SLM.is_onek = 1; % select for new Meadowlark SLM

%%--- Setup SLM Parameters ---%%

if Setup.SLM.is_onek
    Setup.SLM.bit_depth = 12; %For the 512L bit depth is 16, for the small 512 bit depth is 8
else
    Setup.SLM.bit_depth = 8;
end

Setup.SLM.num_boards_found = libpointer('uint32Ptr', 0);
Setup.SLM.constructed_okay = libpointer('int32Ptr', 0);
Setup.SLM.is_nematic_type = 1;
Setup.SLM.RAM_write_enable = 1;
Setup.SLM.use_GPU = 0;
Setup.SLM.max_transients = 10;
Setup.SLM.external_Pulse = 1;
Setup.SLM.timeout_ms = 5000;
% Setup.SLM.lut_file = 'C:\Users\Holography\Desktop\meadowlark\LUT Files\slm6257_at1035_1st_order.lut';
Setup.SLM.lut_file = 'C:\Users\Holography\Desktop\meadowlark\LUT Files\slm6257_at1035_0th_order.lut';
Setup.SLM.reg_lut = libpointer('string');
Setup.SLM.true_frames = 3;
Setup.SLM.pixelmax = 255;%190; changed 4/10/19 improves diffraction efficiency

if Setup.SLM.is_onek
    Setup.SLM.Nx=1024;
    Setup.SLM.Ny=1024;
    Setup.psSLM = 17e-6; % meters, SLM pixel dimensions
else
    Setup.SLM.Nx = 1920;
    Setup.SLM.Ny = 1152;
    Setup.psSLM = 9.2e-6; % meters, SLM pixel dimensions
end

Setup.SLM.wait_For_Trigger = 0; % Set to 1 before initialization as needed.
Setup.SLM.State = 0;


%%--- CGH Computation Parameters ---%%
Setup.CGHMethod = 2;         % Select 1 for superoposition, 2 for GGS, 3 for novocgh, 4 for 2P NovoCGH
Setup.verbose = 1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
Setup.lambda = 1.03e-6;      % meters    Wavelength of the light
Setup.focal_SLM = 0.15;      % meters    focal length of the telescope lens after slm.
Setup.Nx = Setup.SLM.Nx;     % int       Number of pixels in X direction
Setup.Ny = Setup.SLM.Ny;     % int       Number of pixels in Y direction
Setup.useGPU = 1;            % 1 or 0    Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
Setup.maxiter = 50;          % int       Number of iterations (for all methods explored)
Setup.GSoffset = 0.0;        % 11/26 needs to be zero for globalGS. float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms

% % Specify Low and High threshold for threshold-based cost functions
% NormOptions.HighThreshold = 0.5;
% NormOptions.LowThreshold = 0.1;

% Specify Illumination pattern at the SLM, Uniform here, but tunable in general.
Setup.intensity = 1;
Setup.source = sqrt(Setup.intensity)*(1/(Setup.Nx* Setup.Ny))*ones(Setup.Nx, Setup.Ny);


%%--- Etc. Params ---%%

Setup.Datapath = 'Calib_Data';
Setup.Displaypath = 'Calib_Displays';
Setup.Sutterport = 'COM3';