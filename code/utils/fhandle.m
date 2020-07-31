function out = fhandle(in,mode,op_fw,op_bw)

    switch mode
        case 1
            out = op_fw(in);
        case 2
            out = op_bw(in);
    end
    
    out = out(:);
    
end