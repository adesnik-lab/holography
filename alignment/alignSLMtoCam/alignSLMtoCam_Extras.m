%% view fine multi-search images


% for dataUZ4... basler images {plane}{holo}(y,x,z)
for tmp_this_plane=1:numel(dataUZ3)
    for tmp_this_holo=1:numel(dataUZ3{tmp_this_plane})
        tmp_im = max(dataUZ3{tmp_this_plane}{tmp_this_holo},[],3);
%         tmp_im = mean(dataUZ3{tmp_this_plane}{tmp_this_holo},3);
        figure(6767)
        imagesc(tmp_im)
        title(['Plane ' num2str(tmp_this_plane) ', Holo ' num2str(tmp_this_holo)])
        colorbar
        pause
    end
end

%% shoot a holo

i = 3; % plane
j = 1; % holo


target_ref = targListIndiv{i}{j}(1);
expected_z = xyzLoc(3, target_ref);

disp(['Expected z: ' num2str(expected_z) ' um'])

slm.feed(multiHolos{i}(j,:,:))
multi_pwr = size(slmMultiCoordsIndiv{i}{j},2) * pwr;
laser.set_power(multi_pwr)
bas.preview()
laser.set_power(0)

%% single spot from above

idx_holo = 2;
this_holo = slmMultiCoordsIndiv{i}{j}(:,idx_holo)';
[ sholo, ~, ~ ] = function_Make_3D_SHOT_Holos(Setup,this_holo);
slm.feed(sholo)
laser.set_power(pwr)
bas.preview()
laser.set_power(0)

%%
use_power = 10;

slmCoordsTemp = [[0.4 .4 0 1];...
                 [0.8 .2 0 1];...
                 [0.7 .55 0 1]
                 [0.33 .57 0 1]
                 [0.66 .32 0 1]];...
[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos(Setup, slmCoordsTemp);

slm.feed(HoloTemp)
disp('sent SLM')
figure(124)
clf
imagesc(HoloTemp)
colorbar
title('Hologram sent to SLM')

laser.set_power(use_power)
% pause
bas.preview()
laser.set_power(0)

%%
use_power = 2;

slmCoordsTemp = [0.4 .4 0 1];...
[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos(Setup, slmCoordsTemp);

slm.feed(HoloTemp)
disp('sent SLM')
figure(124)
clf
imagesc(HoloTemp)
colorbar
title('Hologram sent to SLM')

laser.set_power(use_power)
% pause
bas.preview()
laser.set_power(0)

