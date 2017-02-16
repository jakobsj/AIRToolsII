function [X,info,restart] = cav(varargin)

[X,info,restart] = sirt('cav',varargin{:});