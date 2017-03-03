function [A,b] = purge_rows(A,b,Nthr)
%PURGE_ROWS  Remove zero/"small" rows of A and corresponding elements of b.
%
% [A,b] = purge_rows(A,b)
% [A,b] = purge_rows(A,b,Nthr)
%
% Identifies zero rows of the coefficient matrix A and removes them.
% If a right-hand side b is present, the corresponding elements of
% b are also removed (b can be a matrix of several right hand sides). A
% needs to be matrix, i.e., cannot be a function handle representing the
% forward operator.
%
% If a positive Nthr is given as the third argument, then all rows in
% which the number of nonzeros is less than or equal to Nthr will be
% removed.
%
% Use this function to 'clean up' a discretized tomography problem.
% Zero rows do not contribute to the reconstruction.
% Rows with few nonzero elements correspond to pixels near the corners of
% the image, whose reconstructions are highly sensitive to noise.

% Per Chr. Hansen, October 1, 2014, DTU Compute.

if nargin<3, Nthr = 0; end

s = sum(A~=0,2);
I = find(s>Nthr);
A = A(I,:);

if nargin>1 && ~isempty(b), b = b(I,:); end