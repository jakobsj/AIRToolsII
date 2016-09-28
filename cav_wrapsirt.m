function [X,info,restart] = cav_wrapsirt(varargin)

[X,info,restart] = sirt('cav',varargin{:});