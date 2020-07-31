function [vecInterp_full_left, vecInterp_full_right] = create_full_vec(alpha,colFull, rowFull,N,d,pixelSize,colAng)

% create full version vectors at farthest depth plane (alpha = 1)

%% resize the full PSF 


% find center of nonzero elements
nz_left = find(colFull~=0);
nz_right = find(rowFull~=0);

center_left = nz_left(round(length(nz_left)/2));
center_right = nz_right(round(length(nz_right)/2));


% d = 4e3; % MASK-TO-SENSOR DISTANCE
% pixelSize = 5.86*2; % PHYSICAL SIZE OF EACH PIXEL MEASUREMENT (IF BINNING, THEN IT WILL BE THE SIZE OF BINNED PIXEL)
% colAng = 25; 
rowAng = colAng;
% N = 128; 
colShifts = -tand(linspace(-colAng,colAng,N))*d/pixelSize;
rowShifts = -tand(linspace(-rowAng,rowAng,N))*d/pixelSize;

center_col = numel(colFull)/2;% nz_left(round(length(nz_left)/2));
center_row = numel(rowFull)/2;% nz_right(round(length(nz_right)/2));
sensorLocations_col = (1:length(colFull)) - center_col;
sensorLocations_row = (1:length(rowFull)) - center_row;

vecInterp_full_left = interp1(sensorLocations_col(:)*alpha,colFull, sensorLocations_col(:),'linear',0);
vecInterp_full_right = interp1(sensorLocations_row(:)*alpha, rowFull, sensorLocations_row(:),'linear',0); 



% % resize PSFs according to the center
% vecInterp_left = imresize(colFull(nz_left), [round(1/alpha*length(nz_left)),1]);
% vecInterp_right = imresize(rowFull(nz_right), [round(1/alpha*length(nz_right)),1]);
% 
% interpLen_left = length(vecInterp_left);
% interpLen_right = length(vecInterp_right);
% 
% % put the resized vectors back to the full vectors
% vecInterp_full_left = zeros(size(colFull));
% vecInterp_full_right = zeros(size(rowFull));
% 
% vecInterp_full_left(center_left -floor(interpLen_left/2) + [0:interpLen_left-1]) = vecInterp_left;
% vecInterp_full_right(center_right -floor(interpLen_right/2) + [0:interpLen_right-1]) = vecInterp_right;


end