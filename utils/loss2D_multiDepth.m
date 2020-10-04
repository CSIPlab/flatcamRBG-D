function [J,grad] = loss2D_multiDepth(depth, est,y,opts)
D = opts.D;
est.depth = reshape(depth,opts.Nd,opts.Nd);
vec = @(x) x(:);
sigma = opts.sigma;
lambda = opts.lambda;

[msr,f1,f2] = fast_msr(est,opts);
[g1,g2] = fast_grad(est,opts, f1, f2);

R = msr - y;
J = norm(R(:))^2; % loss function 
grad = reshape(2*est.image(:)'.*(sum(g1.*(R*f2))+sum(f1.*(R*g2))),[opts.N,opts.N]);

% add TVl2 regularization

% --- Weighted TVL2 ---
W1 = exp(-(D*est.depth).^2/sigma);
W2 = exp(-(est.depth*D').^2/sigma);

lambda_d = 5e-6; 
J = J + lambda*(sum(vec(W1.*(D*est.depth).^2)) + sum(vec(W2.*(est.depth*D').^2)))-lambda_d*sum(log(1-est.depth(:))); 
grad = grad + 2*lambda*(D'*(W1.^2.*(D*est.depth)) + (W2.^2.*(est.depth*D'))*D)+lambda_d./(1-est.depth);


figure(111); 
subplot(121); imagesc(est.depth); title('estimated depth'); colormap jet; colorbar
subplot(122); imagesc(grad); title('gradient'); colorbar
drawnow; 
grad = grad(:);



if opts.Nd == 1; grad = sum(grad(:)); end
end