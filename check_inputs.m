function [Afun,b,m,n,K,Knew,kmax,x0,nonneg,boxcon,L,stoprule,taudelta, ...
    lambdainput,s1,M,w,res_dims,ncp_smooth] = check_inputs(A,b,K,x0,options)

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


%% After here, fix options etc.
if nargin < 5
    options = struct;
end

% Nonnegativity. Default false.
if isfield(options,'nonneg')
    nonneg = options.nonneg;
else
    nonneg = false;
end

% Box constraints [0,L].
if isfield(options,'ubound')
    nonneg = true;
    boxcon = true;
    L = options.ubound;
else
    boxcon = false;
    L = nan;
end

% Default stoprules
stoprule = 'NONE';
taudelta = nan;
if isfield(options,'stoprule') && isfield(options.stoprule,'type')
    stoprule = options.stoprule.type;
end
if isfield(options,'stoprule') && isfield(options.stoprule,'taudelta')
    taudelta = options.stoprule.taudelta;
end

lambdainput = nan;
if isfield(options,'lambda')
    lambdainput = options.lambda;
end

s1 = nan;
if isfield(options,'restart') && isfield(options.restart,'s1')
    s1 = options.restart.s1;
end

% If M is given as input
M = nan;
if isfield(options,'restart') && isfield(options.restart,'M')
    M = options.restart.M;
end

% If weights are given as input
w = nan;
if isfield(options,'w')
    w = options.w;
end

res_dims = n;
if isfield(options,'stoprule') && isfield(options.stoprule,'res_dims')
    res_dims = options.stoprule.res_dims;
    if ~strcmp(stoprule,'NCP')
        warning(['options.stoprule.dims given but only used by NCP ',...
            'stoprule, not the currently specified.']);
    end
end

ncp_smooth = 4;
if isfield(options,'stoprule') && isfield(options.stoprule,'ncp_smooth')
    ncp_smooth = options.stoprule.ncp_smooth;
    if ~strcmp(stoprule,'NCP')
        warning(['options.stoprule.dims given but only used by NCP ',...
            'stoprule, not the currently specified.']);
    end
end
