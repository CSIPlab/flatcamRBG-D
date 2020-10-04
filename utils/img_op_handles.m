function out = img_op_handles(in,mode,f,ft)
switch mode
    case 1
        out = f(in);
    case 2
        out = ft(in);
end
end