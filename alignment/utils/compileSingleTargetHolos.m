function hololist = compileSingleTargetHolos(Setup, slmCoords)

parfor i=1:size(slmCoords, 2)
    t=tic;
    fprintf(['Holo ' num2str(i)]);
    subcoordinates = slmCoords(:,i);
    [Hologram,~,~] = function_Make_3D_SHOT_Holos(Setup,subcoordinates');
    hololist(:,:,i) = Hologram;
    fprintf([' took ' num2str(toc(t)) 's\n']);
end