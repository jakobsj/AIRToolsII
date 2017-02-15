function [X,info,restart] = sart_wrapsirt(varargin)

% Note that any s1 input given to SART is ignored as the value 1 is used.
%
% Note that it is not possible to use weights with SART, as this would
% potentially change the largest singular value away from the known value
% of 1, thus prompting need to compute it numerically, which is now a
% feature of SART to avoid this potentially expensive computation.

[X,info,restart] = sirt('sart',varargin{:});