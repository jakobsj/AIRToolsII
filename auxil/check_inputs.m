function [Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, ...
    relaxparinput,rho,res_dims,rkm1,dk,do_waitbar,verbose,damp,THR,...
    Kbegin,Nunflag] = check_inputs(A,b,K,x0,options)
%CHECK_INPUTS  Check inputs and set default values
%
%   [Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, ...
%    relaxparinput,rho,w,res_dims,rkm1,dk,damp,THR,Kbegin,Nunflag] = ...
%    check_inputs(A,b,K,x0,options)
%
% From inputs given by user to an ART, CART or SIRT method check for
% inconsistencies, extract options or set default values.
% 
% Input:
%    A              Forward operator matrix or function.
%    b              Right-hand side vector.
%    K              List of iterations to return.
%    x0             Starting vector.
%    options        Struct with options.
%    
% Output:
%    Afun           Function handle to forward operator matrix or function.
%    b              Right-hand side vector.
%    m              Number of rows in A.
%    n              Number of columns in A.
%    K              List of iterates to return.
%    kmax           Maximal number of iterations to run.
%    x0             Starting vector.
%    lbound         Lower bound (scalar or vector) on iterates.
%    ubound         Upper bound (scalar or vector) on iterates.
%    stoprule       The name of the stopping rule chosen.
%    taudelta       Stopping rule parameter used by DP and ME.
%    relaxparinput  Relaxation parameter specified by user or default.
%    rho            Spectral radius of iteration matrix.
%    res_dims       Dimensions of residual needed for NCP stopping rule.
%    rkm1           Initialized variable to hold residual at previous 
%                   iteration.
%    dk             Initizlized vector for averaging in NCP.
%    do_waitbar     Whether to display progress using waitbar.
%    verbose        Scalar stating frequency of printing progress.
%    damp           Scalar stating level of damping in ART and CART 
%                   methods.
%    THR            Threshold for use in CART.
%    Kbegin         Flagging parameter for CART.
%    Nunflag        Unflagging parameter for CART.
%
% See also: art, cart, sirt

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Check first number and dimensions of inputs given are correct.

% Check that at least 3 inputs are given.
if nargin < 3
    error('Too few input arguments.')
end

% If A is not a function (i.e., A is a matrix or an object), convert A
% to a function.
if isa(A,'function_handle')
    Afun = A;
else
    Afun = @(xx,transp_flag) afun_matrix(xx,transp_flag,A);
end

% Check that the sizes of A and b match.
mn = Afun([],'size');
m = mn(1); n = mn(2);
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match.')
end

% Check that K is specified and initialize X matrix.
if isempty(K)
    error('Max no. of iterations must ALWAYS be specified.')
end

% Elements of K must appear in strictly increasing order.
if ~all(diff(K)>0)
    error('Elements of K must be strictly increasing.')
end
kmax = K(end);

% Default value for x0.
if nargin < 4 || isempty(x0)
    x0 = zeros(n,1);
end

% Check the size of x0.
if size(x0,1) ~= n || size(x0,2) ~= 1
    error('The size of x0 does not match the problem.')
end


% After here, extract from options or use defaults.
if nargin < 5
    options = struct;
end


% Bound constraints lbound and ubound.

% Lower bound(s). Default empty. Can be specified as either scalar or
% vector of same length as x. If vector, check length/orientation.
if isfield(options,'lbound')
    lbound = options.lbound(:);
    if length(lbound) ~= n && length(lbound) ~= 1
        error(['lbound (lower bound) must either be scalar or ',...
            'same length as vector of unknowns x.'])
    end
else
    lbound = [];
end

% Upper bound(s). Default empty. Can be specified as either scalar or
% vector of same length as x. If vector, check length/orientation.
if isfield(options,'ubound')
    ubound = options.ubound(:);
    if length(ubound) ~= n && length(ubound) ~= 1
        error(['ubound (upper bound) must either be scalar ',...
            'or same length as vector of unknowns x.'])
    end
else
    ubound = [];
end


% Relaxation parameter.
relaxparinput = [];
if isfield(options,'relaxpar')
    relaxparinput = options.relaxpar;
end


% Largest singular value.
rho = nan;
if isfield(options,'rho')
     rho = options.rho;
end


% Stopping rules.

% Default stopping rule.
stoprule = 'NONE';
if isfield(options,'stoprule') && isfield(options.stoprule,'type')
    stoprule = options.stoprule.type;
end

% For stopping rules DP and ME, taudelta must be given in options.
taudelta = nan;
if strcmpi(stoprule,'DP') || strcmpi(stoprule,'ME')
    if isfield(options,'stoprule') && ...
       isfield(options.stoprule,'taudelta') && ...
       isscalar(options.stoprule.taudelta) && ...
       isnumeric(options.stoprule.taudelta)
        taudelta = options.stoprule.taudelta;
    else
        error(['For stopping rules DP and ME a scalar ',...
            'taudelta must be given in options.'])
    end
end

% For stopping rule NCP, res_dims must be given in options.
if strcmp(stoprule,'NCP') && ~isfield(options.stoprule,'res_dims')
    error(['options.stoprule.res_dims must be specified for ',...
        'NCP stopping rule.'])
end
res_dims = []; 
if isfield(options,'stoprule') && isfield(options.stoprule,'res_dims')
    res_dims = options.stoprule.res_dims;
    if ~strcmp(stoprule,'NCP')
        warning(['options.stoprule.res_dims given but only used by NCP ',...
            'stopping rule, not the one currently specified.']);
    end
end

% For stopping rule NCP, set filter length from options or to default.
ncp_smooth = 2;
if isfield(options,'stoprule') && isfield(options.stoprule,'ncp_smooth')
    ncp_smooth = options.stoprule.ncp_smooth;
    if ~strcmp(stoprule,'NCP')
        warning(['options.stoprule.ncp_smooth given but only used by NCP ',...
            'stopping rule, not the one currently specified.']);
    end
end

% Initialize stopping rules
rkm1 = nan; % Only used in ME , otherwise nan as placeholder.
dk = nan;   % Only used in NCP, otherwise nan as placeholder.

switch upper(stoprule)
    case 'DP'
        % DP stopping rule: Nothing to do.
        
    case 'ME'
        % ME stopping rule.
        rkm1 = nan(m,1);
        
    case 'NCP'
        % NCP stopping rule.
        dk = inf(ncp_smooth,1);
        
    case 'NONE'
        % No stopping rule: Nothing to do.
        
    otherwise
        error('The chosen stopping rule is not valid.');
end

% ME is only available for SIRT. If specified for ART or CART, throw error.
stack = dbstack;
if strcmpi(stoprule,'ME') && ~strcmp(stack(2).name,'sirt')
    error(['Stopping rule ME is only available for SIRT methods, ',...
        'not ART and CART.'])
end


% Damping in ART and CART methods.
damp = 0;
if isfield(options,'damp')
    damp = options.damp;
    if damp<0, error('Damping must be positive'), end
end

% CART flagging.

THR = 1e-4;
if isfield(options,'THR')
    THR = options.THR;
end

Kbegin = 10;
if isfield(options,'Kbegin')
    Kbegin = options.Kbegin;
end

Nunflag = round(kmax/4);
if isfield(options,'Nunflag')
    Nunflag = options.Nunflag;
end

% Waitbar and printing info during run.

do_waitbar = false;
if isfield(options,'waitbar') && options.waitbar
    do_waitbar = options.waitbar;
end

verbose = 0;
if isfield(options,'verbose') && options.verbose
    verbose = options.verbose;
end