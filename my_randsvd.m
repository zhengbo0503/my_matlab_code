function [A,evals] = my_randsvd(n, kappa, mode)
%MY_RANDSVD Customized gallery('randsvd')
% Generate symmetric matrices with predefined condition number
% and different distribution of singular values.
% Input :
%      n : dimension of the output
%  kappa : pre-defined 2-norm condition number for the output,
%          if kappa is negative, it will return symmetric p.d. matrix.
%		   Default value is 100. 
%   mode : "geo" or "ari", representing geometrically distributed
%          and arithmetically distributed singular values respectively.  
%          Default value is 'geo'.

classname = 'double';
switch nargin
	case 1
		kappa = 100;
		mode = 'Geo';
	case 2
		mode = 'Geo';
end

if isempty(mode), mode = 'Geo'; end
if sign(kappa) == 1, pd = false; else, pd = true; end
kappa = abs(kappa);

switch mode % Set up vector sigma of singular values.
  case {'Geometric','geometric','Geo','geo'} 
    % geometrically distributed singular values
    factor = kappa^(-1/(n-1));
    sigma = factor.^(0:n-1);
  case {'Arithematic', 'arithematic', 'Ari', 'ari'}
    % arithmetically distributed singular values
    sigma = ones(n,1) - (0:n-1)'/(n-1)*(1-1/kappa); sigma = sigma';
  otherwise
    error(message('MATLAB:randsvd:invalidMode'));
end

sigma = cast(sigma,classname);

if pd
    signs = sign(rand(1,n-2));
    sigma(2:n-1) = sigma(2:n-1).* signs;
end
Q = qmult(n,1,classname);
A = Q*diag(sigma)*Q'; 
A = (A + A')/2; % ensure symmetry

% output true eigenvalues/singular values
if nargout == 2
	evals = sigma;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function above rely on the following function

function B = qmult(A,method,classname)
%QMULT Pre-multiply matrix by random orthogonal matrix.
%   QMULT(A) returns Q*A where Q is a random real orthogonal matrix
%   from the Haar distribution of dimension the number of rows in A.
%   Special case: if A is a scalar then QMULT(A) is the same as QMULT(EYE(A)).
%   QMULT(A,METHOD) specifies how the computations are carried out.
%   METHOD = 0 is the default, while METHOD = 1 uses a call to QR,
%   which is much faster for large dimensions, even though it uses more flops.

%   Called by RANDCOLU, RANDCORR, RANDJORTH, RANDSVD.

%   Reference:
%   G. W. Stewart, The efficient generation of random
%   orthogonal matrices with an application to condition estimators,
%   SIAM J. Numer. Anal., 17 (1980), 403-409.
%
%   Nicholas J. Higham
%   Copyright 1984-2020 The MathWorks, Inc.

if nargin < 2 || isempty(method)
    method = 0;
end
%  Handle scalar A.
if isscalar(A)
   n = A;
   A = eye(n,classname);
else
    n = size(A, 1);
end
if isempty(A) % nothing to do.
    B = A;
    return
end
if method == 1
   [Q,R] = qr(randn(n));
   B = Q*diag(sign(diag(R)))*A;
   return
end
d = zeros(n,1,classname);
for k = n-1:-1:1
    % Generate random Householder transformation.
    x = randn(n-k+1,1);
    s = norm(x);
    sgn = mysign(x(1));
    s = sgn*s;
    d(k) = -sgn;
    x(1) = x(1) + s;
    beta = s*x(1);
    % Apply the transformation to A.
    y = x'*A(k:n,:);
    A(k:n,:) = A(k:n,:) - x*(y/beta);
end
% Tidy up signs.
for i=1:n-1
    A(i,:) = d(i)*A(i,:);
end
A(n,:) = A(n,:)*mysign(randn);
B = A;
end

function S = mysign(A)
%MYSIGN True sign function with MYSIGN(0) = 1.

%   Called by various matrices in elmat/private.
%
%   Nicholas J. Higham
%   Copyright 1984-2013 The MathWorks, Inc.
S = sign(A);
S(S==0) = 1;
end


