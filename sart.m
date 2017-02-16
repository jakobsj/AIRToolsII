function [X,info,restart] = sart(varargin)

% Note that any s1 input given to SART is ignored as the value 1 is used.
%
% Note that it is not possible to use weights with SART, as this would
% potentially change the largest singular value away from the known value
% of 1, thus prompting need to compute it numerically, which is now a
% feature of SART to avoid this potentially expensive computation.
%
% Note that the matrix-free version of SART assumes that all entries in the
% system matrix are nonnegative. Using a system matrix with negative
% entries leads to undefined behavior. In the case of providing the matrix
% explicitly, negative entries are allowed and handled, but this is not
% possible without explicit access to system matrix entries.

[X,info,restart] = sirt('sart',varargin{:});