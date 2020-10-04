function [msr,f_col,f_row] = fast_msr(est,opts)

% est object has two variables. est.image (NxNx3) and est.depth (NxN)

% stack all the row angles into a vector
colShifts_stack = repmat(opts.colShifts',1,opts.N);
% compute f for row dimension
intersect_col = opts.sensorLocations_col'*est.depth(:)' + colShifts_stack(:)';
f_col = interp1(opts.sensorLocations_col(:),opts.colFull(:),intersect_col,'linear', 0);
f_col = f_col(opts.crop_start_col+[1:opts.M],:);

% stack all the col angles into a vector
tmp = repmat(opts.rowShifts,opts.N,1);
rowShifts_stack = tmp(:)';
% compute f for col dimension
intersect_row = opts.sensorLocations_row'*est.depth(:)' + rowShifts_stack;
f_row = interp1(opts.sensorLocations_row(:),opts.rowFull(:),intersect_row,'linear', 0);
f_row = f_row(opts.crop_start_row+[1:opts.M],:);

% compute measurements
image_vec = est.image(:); % stack image into a vector
msr = (f_col.*(image_vec'))*f_row';

end