function y = afun_matrix(x,transp_flag,A)
%AFUN_MATRIX Wraps matrix to allow calling through function handle
%
%   y = afun_matrix(x,transp_flag,A)
%
% This function takes a matrix and wraps it into a function, which can be
% used to apply the matrix in a "pseudo" matrix-free way. Dependent on the 
% second input, either A or A^T is multiplied onto the first input vector 
% (forward or backward operation), or the size of A is returned.
%
% Typical usage is given a matrix A to wrap it in an anonymous function
% myfun = @(XX,TT) afun(XX,TT,A);
% after which myfun implicitly can apply A or A^T or return size of A by
%
% y = myfun(x,'notransp');
% z = myfun(y,'transp');
% s = myfun([],'size');
%
% Input:
%   x           Vector on which to apply matrix multiplication from the
%               left by either A or A^T. x must be a column vector with
%               length matching the relevant dimension of A.
%   transp_flag String to indicate whether to apply forward
%               operation/multiplication by A ('notransp'), backward
%               operation/multiplication by A^T ('transp') or return the
%               size of A ('size'). If set to 'size', the first input is
%               ignored.
%   A           The matrix to be wrapped inside the function.
%
% Output:
%   y           If transp_flag is 'notransp', y is A*x. If 'transp', y is 
%               A'*x. If 'size', y will be size(A).
%
% See also: demo_matrixfree.

% Jakob Sauer Jorgensen.
% 2017-03-03 DTU Compute.

switch transp_flag
    case 'size'
        y = size(A);
    case 'notransp'
        y = A*x;
    case 'transp'
        y = A'*x;
    otherwise
        error('transp_flag must be ''size'', ''notransp'' or ''transp''')
end
