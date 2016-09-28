function [X,info,restart] = landweber_wrapsirt(varargin)

[X,info,restart] = sirt('landweber',varargin{:});