function checkAdjoint(f_h,x)
 x = randn(size(x));
 Ax = f_h(x,1);
 y = randn(size(Ax));
 Aty = f_h(y,2);
 
 innerProduct1 = Ax(:)'*y(:);
 innerProduct2 = x(:)'*Aty(:);
 error = abs(innerProduct1-innerProduct2)...
     /max(abs(innerProduct1),abs(innerProduct2));
 assert(error<1e-9,'"At" is not the adjoint of "A".  Check the definitions of these operators.'); 
end 
