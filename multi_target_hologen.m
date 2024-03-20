function holos = multi_target_hologen(xyzLoc, slmCoordsInitial, opts)


finePts = 13; % odd number please
fineRange = 40;
coaseIncludeThreshScalar = 3;
nframesCapture = 10;

iterationsBeforeStop = 1000;
distanceThreshold = 30; % spacing between points, changed from 50 on 7/15/20 bc new cam
size_of_holo = 10; % number of points in holo
cores = 2; % for parallel compute of holos

zdepths = unique(xyzLoc(3,:));
n_planes = numel(zdepths);

clear slmMultiCoords basCoords targ_list

for i=1:n_planes
    z = zdepths(i);
    targ_idx = find(xyzLoc(3,:) == z);
    slmMultiCoords{i} = slmCoordsInitial(:,targ_idx);
    basCoords{i} = xyzLoc(:,targ_idx);
    targ_list{i} = targ_idx;
end

