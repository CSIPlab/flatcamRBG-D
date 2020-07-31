%% presettings
clear;

addpath(genpath('utils'));

% The list of 4 example objects is: single_card, two_cards, cup and hand
load(fullfile('..','dataset','two_cards'));

d = 0.4; % sensor-to-mask distance (cm)
pixelSize = 5.86*2*1e-4; % sensor pitch (cm)

z2a = @(z) (1-d./z);
a2z = @(a) d./(1-a);

N = 128; M = 512; npatterns = N;
measure_row = 512; measure_col = 512;
colAng = 30; rowAng = 30;

vec = @(x) x(:);
nor1 = @(x) x/max(abs(x(:)));
mean_sub = @(x) (x - ones(M,1)*mean(x,1) - mean(x,2)*ones(1,M) + mean(mean(x)));

reshp = @(x) reshape(x,N,N);

eval dataProcessing_flashlight

[colFull_inf, rowFull_inf] = create_full_vec(z2a(42),colFull,rowFull,N,d,pixelSize,colAng); % create full version vectors when the PSF is at infinite far

diff = [1,-1]';
colFull_grad = [0; conv(colFull_inf, diff,'valid')];
rowFull_grad = [0; conv(rowFull_inf, diff,'valid')];

nz_col = find(colFull_inf~=0);
nz_row = find(rowFull_inf~=0);
center_col = numel(colFull)/2;
center_row = numel(rowFull)/2;

sensorLocations_col = (1:length(colFull_inf)) - center_col;
sensorLocations_row = (1:length(rowFull_inf)) - center_row;


colShifts = -tand(linspace(-colAng,colAng,N))*d/pixelSize;
rowShifts = -tand(linspace(-rowAng,rowAng,N))*d/pixelSize;

opts = [];
opts.d = d;
opts.pixelSize = pixelSize;
opts.M = M; opts.N = N;
opts.crop_start_col = 512;
opts.crop_start_row = 512;
opts.colFull = colFull_inf;
opts.rowFull = rowFull_inf;
opts.colFull_grad = colFull_grad;
opts.rowFull_grad = rowFull_grad;
opts.sensorLocations_col = sensorLocations_col;
opts.sensorLocations_row = sensorLocations_row;
opts.colShifts = colShifts;
opts.rowShifts = rowShifts;
opts.D = convmtx([1,-1],N-1);

%% Parameters

% LSQR parameters
damp = 0.4; % l2 regularization magnitude
atol   = 1.0e-6;
btol   = 1.0e-6;
conlim = 1.0e+8;
itnlim = 50;
show   = 0;

% minFunc parameters
options = [];
options.display = 10; % 'none';
options.maxFunEvals = 50;
options.maxIter = 20;
options.Method = 'lbfgs';
options.progTol = 1e-8;
options.useMex = 1;

niter = 5;

res_images = zeros(N,N,3,1);
res_depths = zeros(N,N,1);
res_images_gdy_raw = zeros(N,N,3,1);
res_images_gdy = zeros(N,N,3,1);
res_depths_gdy = zeros(N,N,1);

%% select scenes
jj = 1;
Y = mean(realobj_table(:,:,:,:,jj),4)-background_image_realobj; % set target measurements

est.image = zeros(N);
est.depth = ones(N)*z2a(20); % initialize estimated image and depth

%% create stack of A
% create a depth list
ndepth = 10;
depth_list = linspace(z2a(15),z2a(40),ndepth);

reshp_gdy = @(x) reshape(x,N,N,ndepth);

% create stack of transfer matrices
A_left_dl = zeros(M,N,ndepth);
A_right_dl = zeros(M,N,ndepth);

for dd = 1:ndepth
    [A_left_1, A_right_1] = eval_PSF(depth_list(dd),M,N,opts);
    A_left_dl(:,:,dd) = A_left_1;
    A_right_dl(:,:,dd) = A_right_1;
end

%% depth sweep
image_depthSwipe = zeros(N,N,3,ndepth);
loss_depthSwipe = zeros(1,ndepth);

bestAlpha = 1;
bestLoss = inf;
for dd = 1:ndepth
    [A_left_1, A_right_1] = eval_PSF(depth_list(dd),M,N,opts);
    A_left_dl(:,:,dd) = A_left_1;
    A_right_dl(:,:,dd) = A_right_1;
    
    Y = mean(realobj_table(:,:,:,:,jj),4)-background_image_realobj; % set target measurements
    f_img_op_plane = @(x,mode) image_operator(x,mode,A_left_1,A_right_1);
    
    xh1 = reshp(lsqrSOL(M^2,N^2,f_img_op_plane,vec(mean_sub(Y(:,:,1))),damp,atol,btol,conlim,itnlim,show));
    xh2 = reshp(lsqrSOL(M^2,N^2,f_img_op_plane,vec(mean_sub(Y(:,:,2))),damp,atol,btol,conlim,itnlim,show));
    xh3 = reshp(lsqrSOL(M^2,N^2,f_img_op_plane,vec(mean_sub(Y(:,:,3))),damp,atol,btol,conlim,itnlim,show));
    
    xh_RGB = cat(3,nor1((xh1)), nor1((xh2)), nor1((xh3)) );
    image_depthSwipe(:,:,:,dd) = xh_RGB;
    yy = cat(3,mean_sub(Y(:,:,1)),mean_sub(Y(:,:,2)),mean_sub(Y(:,:,3)));
    loss = norm(vec(cat(3,A_left_1*xh1*A_right_1',A_left_1*xh2*A_right_1',A_left_1*xh3*A_right_1') - (yy))) / norm(vec((yy)));
    loss_depthSwipe(dd) = loss;
    
    figure(11);
    imagesc(xh_RGB); title(sprintf('scene %d, depth %0.5gcm, loss ratio %0.5g', jj,a2z(depth_list(dd)),loss));
    drawnow
    
    if loss < bestLoss
        bestLoss = loss;
        bestAlpha = depth_list(dd);
    end
end

%% continuous depth

opts.sigma = inf;
opts.lambda = 50;

% estimate image using current depth
fprintf('Initialization with depth %3.4g \n',a2z(bestAlpha));
est.image = xh2;
est.depth = ones(N)*bestAlpha;
% est.depth = ones(N)*z2a(5);
[~,f1,f2] = fast_msr(est,opts);

y = Y(:,:,2);

f_img_op_allAng = @(x,mode) image_operator_allAng(mode,x,f1,f2,opts);
checkAdjoint(f_img_op_allAng,randn(N^2,1));

x = lsqrSOL(M^2,N^2,f_img_op_allAng,vec(mean_sub(y)),damp,atol,btol,conlim,itnlim,show);
x(x<0) = 0;
x(x<median(x)) = 0;
est.image = reshape(x,N,N);

Nd = N;
opts.Nd = Nd;

error_table = [];

for iter = 1:niter
    % estimate depth
    depth_handle = @(depth) loss2D_multiDepth(depth, est,mean_sub(y),opts);
    est.depth = reshape(minFunc(depth_handle,vec(est.depth),options),Nd,Nd);
    
    % estimate image using current depth
    [~,f1,f2] = fast_msr(est,opts);
    f_img_op_allAng = @(x,mode) image_operator_allAng(mode,x,f1,f2,opts);
    checkAdjoint(f_img_op_allAng,randn(N^2,1));
    x = lsqrSOL(M^2,N^2,f_img_op_allAng,vec(mean_sub(y)),damp,atol,btol,conlim,itnlim,show);
    
    
    figure(22);
    subplot(131); imagesc(est.image); title(sprintf('iter %d', iter)); colormap gray; axis image 
    subplot(132); imagesc(a2z(est.depth)); title('depth (cm)'); colormap jet; colorbar; axis image
    if opts.Nd > 1
        subplot(133); depthTo3D(est.image,a2z(est.depth),colAng,45,45);
    end
    drawnow
    
end

%%
% compute colored images using current depth estimate
[msr,f1,f2] = fast_msr(est,opts);
f_img_op_allAng = @(x,mode) image_operator_allAng(mode,x,f1,f2,opts);
checkAdjoint(f_img_op_allAng,randn(N^2,1));

% estimate colored image
[x1,~,~,J1] = lsqrSOL(M^2,N^2,f_img_op_allAng,vec(mean_sub(Y(:,:,1))),damp,atol,btol,conlim,itnlim,show);
[x2,~,~,J2] = lsqrSOL(M^2,N^2,f_img_op_allAng,vec(mean_sub(Y(:,:,2))),damp,atol,btol,conlim,itnlim,show);
[x3,~,~,J3] = lsqrSOL(M^2,N^2,f_img_op_allAng,vec(mean_sub(Y(:,:,3))),damp,atol,btol,conlim,itnlim,show);

x_color = cat(3,nor1(reshp(x1)),nor1(reshp(x2)),nor1(reshp(x3)));
depth_cm = a2z(est.depth);
figure(33);
subplot(121); imagesc(x_color); title('image from continuous'); axis image off
subplot(122); imagesc(depth_cm); colormap gray; colorbar; title('depth map (cm) from continuous'); axis image off

% results
res_images(:,:,:,jj) = x_color;
res_depths(:,:,jj) = depth_cm;

%% plot and save

figure(24);
imagesc(post_process(x_color));  axis image off
export_fig('image_color');

depth_cm(x_color(:,:,1) < 0.15) = nan;
figure(25);
imagesc(depth_cm,[10,40]); axis image off; colormap jet;
set(gcf,'color','w');
cH = colorbar;
set(cH,'fontsize',25);
export_fig('depth_3D');

save('continuous_depth_data.mat',...
    'res_images','res_depths','image_depthSwipe',...
    'x_color','depth_cm');