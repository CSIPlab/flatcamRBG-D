function  cnorms = colnorms(A,m,n,ntrials,nfigure)

% Test:   cnorms = colnorms(A,m,n,5,1);
% Normal: cnorms = colnorms(A,m,n);   is equivalent to
%         cnorms = colnorms(A,m,n,5,0);
%
% colnorms estimates the 2-norms of the columns of the m x n matrix A
% using several products A'*v, where v is a random vector of +/-1s.
% 
% INPUT:
%   A        A dense or sparse matrix, or an operator.
%            If A is an operator, the statement y = A(mode,v);
%            must return y = A*v  when mode = 1 and v is an n-vector,
%                        y = A'*v when mode = 2 and y is an m-vector.
%            (y will have dimension m or n respectively.)
%            
%   ntrials  Number of iterations the algorithm takes (default 5).
%
%   nfigure  =0: Don't compare with true norms (default 0).
%            >0: Compare with true norms and plot histogram in that figure.
%
% OUTPUT:
%   cnorms:  The estimates of the 2-norms
%
% SOURCE:
%   T.-Y. Chen and J. Demmel. Balancing sparse matrices for computing
%   eigenvalues.  Linear Algebra and Its Applics., 309, 261-287, 2000.
%
% 22 Feb 2007: normEst.m written by Kaustuv:
%              Kaustuv (k a u s t u v a t g m a i l d o t c o m)
%              iCME, Stanford University.
%
% 23 Jan 2008: colnorms.m derived from normEst.m
%              Converted from rownorms to colnorms.
%              State of randn preserved.
%              A may be a matrix or an operator.
%              To treat a matrix A as an operator Aop, do this:
%                 Aop    = @(mode,v) Aprod(mode,v,A);
%                 [m,n]  = size(A);
%                 cnorms = colnorms(Aop,m,n,5,0);
%
% 15 Mar 2009: No need to store lots of vectors.
%
% Maintained by Michael Saunders (saunders@stanford.edu).


if nargin<4 || isempty(ntrials) || ntrials<1
   ntrials = 5;
end

if nargin<5 || isempty(nfigure) || nfigure<1
   nfigure = 0;
end

oper  = isa(A,'function_handle');
state = randn('state');     % Save current state
rand('state',0)

cnorms = zeros(n,1);

for trial = 1:ntrials
  v = ones(m,1);
  v(find(randn(m,1) < 0)) = -1;
  if oper
    y = A(2,v);
  else
    y = A'*v;
  end
  cnorms = cnorms + y.^2;
end

cnorms = sqrt(cnorms) / sqrt(ntrials);

% If requested, compare against true norms.
if nfigure > 0
  truenorms = zeros(n,1);
  v         = zeros(n,1);

  for j = 1:n
    if oper
      v(j) = 1;
      y    = A(1,v);   % y = A*ej
      v(j) = 0;
    else
      y    = A(:,j);
    end
    truenorms(j) = norm(y);
  end

  ratio = cnorms ./ truenorms;
  figure(nfigure)
  hist(ratio)
end

randn('state',state);    % Restore state
return
