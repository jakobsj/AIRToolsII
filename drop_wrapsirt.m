function [X,info,restart] = drop_wrapsirt(varargin)

[X,info,restart] = sirt('drop',varargin{:});