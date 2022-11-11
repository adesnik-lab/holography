function holos = generateRandomHolos(xRange, yRange, zRange, npts)

slmRange = [xRange; yRange; zRange];

% reset random number generator
rng('shuffle');

% generate random points and scale across range
pts = rand([3,npts]);
rngScaler = slmRange(:,2) - slmRange(:,1);
pts = rngScaler.*pts;
pts = pts + slmRange(:,1);

% slmCoords are ([X,Y,Z,P], npts)
holos = ones([4,npts]);
holos(1:3,:) = pts;