function out = post_process(in)

% get channel
[m,n,ncha] = size(in);

out = in;

% thresholding
out(out<0) = 0;
out(out<0) = 0;
out(out<median(out)) = 0;

% out = histeq(out);
% % out = out.*(in(:,:,1)>0.12) + in.*(in(:,:,1)<=0.12);

end