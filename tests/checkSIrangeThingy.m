[X,Y] = meshgrid([0:10:510],[0:10:510]);
test_pts = [X(:), Y(:)];
test_pts(:,3)=30;
output_test = function_Eval3DCoC(CoC.SItoSLM, test_pts);
de_list = computeDEfromList(test_pts', num2cell(1:size(test_pts,1)'), ones(1,size(test_pts,1))');

in_slm_bounds = inpolygon(output_test(:,1),output_test(:,2),[0.05,.95],[0.05,.95]);

figure(666)
scatter(test_pts(:,1), test_pts(:,2),[],'r')
hold on;
scatter(test_pts(in_slm_bounds,1), test_pts(in_slm_bounds,2), [], de_list(in_slm_bounds))
colorbar
