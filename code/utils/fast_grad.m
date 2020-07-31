function [g_col,g_row] = fast_grad(est,opts, f1, f2)

% est object has two variables. est.image (NxNx3) and est.depth (NxN)

% stack all the row angles into a vector
colShifts_stack = repmat(opts.colShifts',1,opts.N);

% compute gradients
D = -convmtx([1,-1],opts.M-1);
D(end+1,end) = 0; 
% 
% D = -convmtx([-1,2,-1]/2,opts.M-1);
% D = D(:,2:end);
% D(1,1:2) = 0; 
% D(end+1,end-1:end) = 0; 

g_col = opts.sensorLocations_col(opts.crop_start_col+[1:opts.M])'.*(D*f1); 
g_row = opts.sensorLocations_row(opts.crop_start_row+[1:opts.M])'.*(D*f2); 

% % compute f for row dimension
% intersect_col = opts.sensorLocations_col'*est.depth(:)' + colShifts_stack(:)';
% % g_col = interp1(opts.sensorLocations_col',opts.colFull_grad(:).*opts.sensorLocations_col(:),intersect_col,'linear', 0);
% g_col = opts.sensorLocations_col(:).*(interp1(opts.sensorLocations_col',opts.colFull_grad(:),intersect_col,'linear', 0));
% g_col = g_col(opts.crop_start_col+[1:opts.M],:);
% 
% % stack all the col angles into a vector
% tmp = repmat(opts.rowShifts,opts.N,1);
% rowShifts_stack = tmp(:)';
% % compute f for col dimension
% intersect_row = opts.sensorLocations_row'*est.depth(:)' + rowShifts_stack;
% % g_row = interp1(opts.sensorLocations_row(:),opts.rowFull_grad(:).*opts.sensorLocations_row(:),intersect_row,'linear', 0);
% g_row = opts.sensorLocations_row(:).*interp1(opts.sensorLocations_row(:),opts.rowFull_grad(:),intersect_row,'linear', 0);
% g_row = g_row(opts.crop_start_row+[1:opts.M],:);


end