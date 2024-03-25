% these are for the sutter
mvShortWait = 0.2;
mvLongWait = 2;

% sockets
laser_socket = 42134;
si_socket = 42135;

% this is laser power in mW, but is dynamically adjusted in the script, so don't include for now
% pwr = 5;

% camera background frames
nBackgroundFrames = 10;

% slm XYZ range
slmXrange = [0.06 0.96];
slmYrange = [0.06 0.94];
slmZrange = [-0.022 0.02];

% initial number of SLM single target points
npts = 250;

% ScanImage planes calibration
nOptotunePlanes = 15;
sutterPlanes    = -25:5:150;
nframesCapture  = 5; % camera average
gridpts = 25;
SIthresholdMod = 1.5; % for exclusions
optotuneMinDepth = -10; % for exclusions
optotuneMaxDepth = 130; % for exclusions
