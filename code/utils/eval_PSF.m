function [A_simu_left, A_simu_right,grad_left,grad_right] = eval_PSF(alpha,M,N,opts)

intersect_col = alpha*repmat(opts.sensorLocations_col',[1,N]) + opts.colShifts; % compute the intersection locations on the mask.
intersect_row = alpha*repmat(opts.sensorLocations_row',[1,N]) + opts.rowShifts; 

A_simu_left = interp1(opts.sensorLocations_col',opts.colFull,intersect_col,'linear',0);
A_simu_right = interp1(opts.sensorLocations_row',opts.rowFull,intersect_row,'linear',0);

% cropping
A_simu_left = A_simu_left(opts.crop_start_col+[1:M],:);
A_simu_right = A_simu_right(opts.crop_start_row+[1:M],:);

% gradients
grad_left = interp1(opts.sensorLocations_col',opts.colFull_grad,intersect_col,'linear',0).*opts.sensorLocations_col';
grad_right = interp1(opts.sensorLocations_row',opts.rowFull_grad,intersect_row,'linear',0).*opts.sensorLocations_row';

grad_left = grad_left(opts.crop_start_col+[1:M],:);
grad_right = grad_right(opts.crop_start_row+[1:M],:);

end