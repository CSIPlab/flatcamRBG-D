function x = lsqrtest( m, n, damp )
%        x = lsqrtest( m, n, damp );
%
% If m = n and damp = 0, this sets up a system Ax = b
% and calls lsqr.m to solve it.  Otherwise, the usual
% least-squares or damped least-squares problem is solved.
%
% lsqraprod.m defines the m x n matrix A.

% 11 Apr 1996: First version for distribution with lsqr.m.
%              Michael Saunders, Dept of EESOR, Stanford University.
% 07 Aug 2002: LSQR's output parameter rnorm changed to r1norm, r2norm.
%              Michael Saunders, Systems Optimization Laboratory,
%              Dept of MS&E, Stanford University.
%-----------------------------------------------------------------------
  
xtrue  = (n : -1: 1)';
iw     = 0;
rw     = 0;
b      = lsqraprod( 1, m, n, xtrue, iw, rw );

atol   = 1.0e-6;
btol   = 1.0e-6;
conlim = 1.0e+10;
itnlim = 10*n;
show   = 1;

[ x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var ]...
  = lsqr( m, n, 'lsqraprod', iw, rw, b, ...
          damp, atol, btol, conlim, itnlim, show );

disp(' ');   j1 = min(n,5);   j2 = max(n-4,1);
disp('First elements of x:');  disp(x(1:j1)');
disp('Last  elements of x:');  disp(x(j2:n)');

y    = lsqraprod( 1, m, n, x, iw, rw );
r    = b - y;
r1   = norm(r);
r2   = norm([r; (-damp*x)]);
disp(' ')
str1 = sprintf( 'r1norm, r2norm (est.)  %10.3e %10.3e', r1norm, r2norm );
str2 = sprintf( 'r1norm, r2norm (true)  %10.3e %10.3e', r1    , r2     );
disp(str1)
disp(str2)

%===================
% End of lsqrtest.m
%===================
