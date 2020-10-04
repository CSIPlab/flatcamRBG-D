function out = image_operator_allAng(mode,in,f1,f2,varargin)
if length(varargin); opts = varargin{1}; end

switch mode
    case 1
        % forward operator
        out = f1.*(in')*f2';
    case 2
        % adjoint operator
        in = reshape(in,[opts.M,opts.M]);
        out = sum((f1'*in)'.*f2);        
end
out = out(:);
end