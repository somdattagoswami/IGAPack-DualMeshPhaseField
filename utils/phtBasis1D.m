function [N, dN]=phtBasis1D(u_hat, Ce, p)

%calculate the shape function and first derivatives
%INPUT: e - element index
%       u_hat - evaluation point in u-direction reference coordinates (from -1 to 1)
%       Ce - local Bezier extraction operator
%       p - polynomial degree
%OUTPUT: N - value of the p+1 shape functions at the evaluation point
%        dN - derivatives of the shape functions
% NOTE: The use of this function incurs a performance penalty and should be
% deprecated
% ------------------------------------------------------------------

%     evaluate 1d shape functions and derivatives
[B, dBxi] = bernstein_basis(u_hat,p);

%form B-spline basis functions using the Bezier extraction operator
N = B*Ce';
dN = dBxi*Ce';
