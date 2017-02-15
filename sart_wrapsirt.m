function [X,info,restart] = sart_wrapsirt(varargin)

% Note that any s1 input given to SART is ignored as the value 1 is used.

[X,info,restart] = sirt('sart',varargin{:});