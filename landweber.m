function [X,info,restart] = landweber(varargin)

[X,info,restart] = sirt('landweber',varargin{:});