function [X,info,restart] = cimmino_wrapsirt(varargin)

[X,info,restart] = sirt('cimmino',varargin{:});