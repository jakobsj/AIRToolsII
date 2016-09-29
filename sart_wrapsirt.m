function [X,info,restart] = sart_wrapsirt(varargin)

[X,info,restart] = sirt('sart',varargin{:});