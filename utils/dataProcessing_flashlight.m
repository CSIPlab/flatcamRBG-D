%% data processing modify
extend_size = 512;
vec_extend = @(x) [zeros(extend_size,1);x;zeros(extend_size,1)];
vec_shrink = @(x) x(extend_size+[1:measure_row]);
fsh_process = @(x) mean_sub(x-background_image_flsh);
normal = @(x) x./sqrt(sum(x.^2)); % vector norm to 1

% extendedSensor_locations = -(extend_size+M)/2:(extend_size+M)/2,

% M = measure_row;
% N = npatterns;

% Full column vector

% colFull = zeros(measure_row*2,1);
k = 5; 
[U,S,V] = svds(fsh_process(fsh_table(:,:,k)),1); % create the first column of left matrix
colFirst = (normal((U*sqrt(S))));

[U,S,V] = svds(fsh_process(fsh_table(:,:,k)),1); % create the first column of right matrix
rowFirst = (normal((sqrt(S)*V)));

colFull = vec_extend(colFirst);
rowFull = vec_extend(rowFirst);

colVecShift_stack = zeros(numel(colFull),9); 
rowVecShift_stack = zeros(numel(rowFull),9);
colVecShift_stack(:,k) = colFull; 
rowVecShift_stack(:,k) = rowFull; 

for ii = [1:k-1 k+1:9]  % merge the selected flashlight into full columns
    % col and row responses
    [U,S,V] = svds(fsh_process(fsh_table(:,:,ii)),1);
    colVec = (normal(U*sqrt(S)));
    rowVec = (normal(sqrt(S)*V));
    % correlation 
    corrMap_col = normxcorr2(colVec,colFirst);
    corrMap_row = normxcorr2(rowVec,rowFirst);
    if abs(min(corrMap_col)) > max(corrMap_col) 
        colVec = -colVec; % SVD gives different positive/negative arrangement for each flashlight
    end
    if abs(min(corrMap_row)) > max(corrMap_row)
        rowVec = -rowVec;
    end
    corrMap_col = normxcorr2(colVec,colFirst);
    corrMap_row = normxcorr2(rowVec,rowFirst);
    
    % compute shift distances
    [~,maxInd_col] = max(abs(corrMap_col));
    [~,maxInd_row] = max(abs(corrMap_row));
    
    shift_dist_col = ceil(length(corrMap_col)/2) - maxInd_col;
    shift_dist_row = ceil(length(corrMap_row)/2) - maxInd_row;
    
%     fprintf('%d-th PSF, shift distance (%d,%d) \n',ii, shift_dist_col,shift_dist_row);
    
    % merge new vectors into full vectors
    colVec = vec_extend(colVec);
    rowVec = vec_extend(rowVec);
    
    colVec_shift = circshift(colVec,-shift_dist_col); % shift back vector
    idx_merge = find((colFull==0).*(colVec_shift~=0));
    idx_mean = find((colFull~=0).*(colVec_shift~=0));
    
    colFull(idx_merge) = colVec_shift(idx_merge);
    colFull(idx_mean) = (colFull(idx_mean) + colVec_shift(idx_mean))/2;

    
    
    rowVec_shift = circshift(rowVec,-shift_dist_row); % shift back vector
    idx_merge = find((rowFull==0).*(rowVec_shift~=0));
    idx_mean = find((rowFull~=0).*(rowVec_shift~=0));
    
    rowFull(idx_merge) = rowVec_shift(idx_merge);
    rowFull(idx_mean) = (rowFull(idx_mean) + rowVec_shift(idx_mean))/2;
    
    
    colVecShift_stack(:,ii) = colVec_shift; 
    rowVecShift_stack(:,ii) = rowVec_shift; 

end

% create full sized column and row vectors of the psf
colFull = sum(colVecShift_stack,2)./sum(double(abs(colVecShift_stack)>0),2); 
rowFull = sum(rowVecShift_stack,2)./sum(abs(rowVecShift_stack)>0,2); 

colFull(isnan(colFull)) = 0; 
rowFull(isnan(rowFull)) = 0; 

% shifting first columns to simulate
A_simu_left = zeros(length(colFull),npatterns);
A_simu_right = zeros(length(rowFull),npatterns);

colShifts = -tand(linspace(-colAng,colAng,N))*d/pixelSize;
rowShifts = -tand(linspace(-rowAng,rowAng,N))*d/pixelSize;

for jj = 1:npatterns
    A_simu_left(:,jj) = imtranslate(colFull,[0,colShifts(jj)]);
    A_simu_right(:,jj) = imtranslate(rowFull,[0,rowShifts(jj)]);
end

% cropping
A_simu_left = A_simu_left(extend_size+[1:measure_row],:);
A_simu_right = A_simu_right(extend_size+[1:measure_row],:);
