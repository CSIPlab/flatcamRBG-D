function depthTo3D(est_image,est_depth,theta,view1,view2)
%%
% est_image:NxNx3
% est_depth:NxN image, in cm
N = size(est_image,1);

[x_theta,y_theta] = meshgrid(linspace(-theta,theta,N));

X = est_depth.*tand(x_theta);
Y = -est_depth.*tand(y_theta);
Z = est_depth;

% find negative positions
if size(est_image,3) > 1
    vec = @(z) z(:)
    neg_pos = find(rgb2gray(abs(est_image))<median(vec(rgb2gray(abs(est_image)))));
else
    neg_pos = [];% find((est_image)<0.13);
end 
X_nng = X(:);
X_nng(neg_pos) = [];
Y_nng = Y(:);
Y_nng(neg_pos) = [];
Z_nng = Z(:);
Z_nng(neg_pos) = [];
img_nng = reshape(est_image,[],size(est_image,3));
img_nng(neg_pos,:) = [];

scatter3(X_nng,Z_nng, Y_nng, 12,img_nng,'filled');
axis([min(X(:)) max(X(:)) 10 30 min(Y(:)) max(Y(:)) ])
% scatter3(X(:),Z(:),Y(:),16,reshape(est_image,[],3));
view(view1,view2);
drawnow; 
end