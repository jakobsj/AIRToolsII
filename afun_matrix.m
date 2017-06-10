function y = afun_matrix(x,transp_flag,A)
%AFUN_MATRIX Wraps matrix to allow calling through function handle
%
%   y = afun_matrix(x,transp_flag,A)
%
% This function takes a matrix and wraps it into a function, which can be
% used to apply the matrix in a "pseudo" matrix-free way. Dependening on
% the  second input, either A or A' is multiplied onto the first input
% vector, or the size of A is returned.
%
% Given a matrix A, the typical use of this function is  to wrap A in an
% anonymous function
%    myfun = @(XX,TT) afun(XX,TT,A);
% after which myfun implicitly can apply A or A^T or return the size of A:
%    y = myfun(x,'notransp');
%    z = myfun(y,'transp');
%    s = myfun([],'size');
%
% Input:
%   x           Vector on which to apply matrix multiplication from the
%               left by either A or A'; x must be a column vector with
%               length matching the relevant dimension of A.
%   transp_flag String to indicate whether to apply multiplication by A
%               ('notransp') or A' ('transp'), or return the size of A
%               ('size'). If set to 'size', the first input is ignored.
%   A           The matrix to be wrapped inside the function.
%
% Output:
%   y           If transp_flag is 'notransp', y is A*x. If 'transp', y is 
%               A'*x. If 'size', y will be size(A).
%
% See also: demo_matrixfree.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, DTU Compute, 2010-2017.

% This file is part of the AIR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package. 
% 
% Copyright 2017 Per Christian Hansen & Jakob Sauer Jorgensen, DTU Compute

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
