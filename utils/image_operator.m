function out = image_operator(in, mode,A_left,A_right)

[M,N] = size(A_left);

switch mode
    case 1
        out = A_left*reshape(in,N,N)*A_right';
        
    case 2
        out = A_left'*reshape(in,M,M)*A_right;
end

out = out(:);

end