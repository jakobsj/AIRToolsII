function [Afun,b,m,n,Knew,kmax,x0] = check_inputs(A,b,K,x0,options)

% Add check of options input, including stopping criteria ones, such as
% taudelta

% Probably change landweber etc inputs to collected all in varargin. This
% allows a variable number of inputs (eg with or without options) to be
% passed on to check_inputs and then just return back to landweber Afun etc
% so that not also A is available. Then also rename check_inputs to
% parse_inputs, since inputs are not really just checked put changed into
% the ones to be used.


% Check that at least 3 inputs are given.
if nargin < 3
    error('Too few input arguments.')
end

% If A is not a function (i.e., A is a matrix or an object), convert A
% to a function.
if isa(A,'function_handle')
    Afun = A;
else
    Afun = @(xx,transp_flag) Afun_matrix(xx,transp_flag,A);
end

% Check that the sizes of A and b match.
mn = Afun([],'size');
m = mn(1); n = mn(2);
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match.')
end

% Check that K is specified and initialize X matrix..
if isempty(K)
    error('Max no. of iterations must ALWAYS be specified.')
end
Knew = sort(K);
kmax = Knew(end);

% Default value for x0.
if nargin < 4 || isempty(x0)
    x0 = zeros(n,1);
end

% Check the size of x0
if size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem.')
end