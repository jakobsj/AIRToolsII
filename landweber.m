function [X,info,ext_info] = landweber(varargin)

[X,info,ext_info] = sirt('landweber',varargin{:});