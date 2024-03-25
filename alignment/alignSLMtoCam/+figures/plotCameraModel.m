function plotCameraModel(fits)

cam2XYZ = fits.cam2XYZ;
obsZ = fits.obsZ;
OptZToCam = fits.OptZToCam;

model_eval = polyvaln(OptZToCam, cam2XYZ(1:3,:)');

newfig('Camera Model Fits')
f = gcf();
f.Units = 'normalized';
f.Position = [0.4 0.4 0.4 0.4];

subplot(2,2,1)
scatter3(cam2XYZ(1,:), cam2XYZ(2,:), cam2XYZ(3,:), [], obsZ, 'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Measured Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,2)
scatter3(cam2XYZ(1,:), cam2XYZ(2,:), cam2XYZ(3,:), [], model_eval, 'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Estimated Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,3)
scatter3(cam2XYZ(1,:), cam2XYZ(2,:), cam2XYZ(3,:), [], model_eval-obsZ, 'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Error (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,4)
c = sqrt((model_eval-obsZ).^2);
scatter3(cam2XYZ(1,:), cam2XYZ(2,:), cam2XYZ(3,:), [], c, 'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Error RMS (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square