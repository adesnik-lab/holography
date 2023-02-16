slmVals = [0:0.01:1; 0:0.01:1];
slmVals(3,:) = -.03;

% approach 1, long way
estCam = function_Eval3DCoC(CoC.SLMtoCam, slmVals');
SIvals = function_Eval3DCoC(CoC.CamToSI, estCam);
estOpto = polyvaln(CoC.camToOpto,estCam);
DEvals = DEfromSLMCoords(slmVals');

% approach 2, aka shortcut
% SIvals = function_Eval3DCoC(CoC.SLMtoSI, slmVals);

figure(420)
clf


% subplot(1,2,1)
scatter(slmVals(1,:), SIvals(:,1), [], DEvals, 'filled')
hold on
scatter(slmVals(2,:), SIvals(:,2), [], DEvals, 'filled')
xlabel('SLM')
ylabel('SI')
legend('X', 'Y')
caxis([0 1])
colorbar

% subplot(1,2,2)
% plot(slmVals(3,:), estOpto, 'LineWidth', 2)


disp('---')
disp(['zero-oder x min: ' num2str(SIvals(slmVals(1,:) == 0.4, 1))])
disp(['zero-oder x max: ' num2str(SIvals(slmVals(1,:) == 0.6, 1))])
disp(['zero-oder y min: ' num2str(SIvals(slmVals(2,:) == 0.45, 2))])
disp(['zero-oder y max: ' num2str(SIvals(slmVals(2,:) == 0.64, 2))])
disp(' ')

line([0.4 0], [SIvals(slmVals(1,:) == 0.4,1) SIvals(slmVals(1,:) == 0.4,1)], 'Color', 'k', 'LineWidth',2)
line([0.6 0], [SIvals(slmVals(1,:) == 0.6,1) SIvals(slmVals(1,:) == 0.6,1)], 'Color', 'k', 'LineWidth',2)

line([0.45 0], [SIvals(slmVals(2,:) == 0.45,2) SIvals(slmVals(2,:) == 0.45,1)], 'Color', 'g', 'LineWidth',2)
line([0.64 0], [SIvals(slmVals(2,:) == 0.64,2) SIvals(slmVals(2,:) == 0.64,1)], 'Color', 'g', 'LineWidth',2)
