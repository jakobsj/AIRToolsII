function [Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, ...
    relaxparinput,s1,w,res_dims,ncp_smooth,damp,THR,Kbegin,Nunflag] = ...
    check_inputs(A,b,K,x0,options)

% PCH: removed "Knew" from output.

% PCH: we ought to check if the user inputs options.stoprule.taudelta
%      is case of DP and ME.

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
% PCH: added this check that the elements of the user's K appear in
%      strictly increasing order:
if ~all(diff(K)>0)
    error('Elements of K must be strictly increasing.')
end
% Knew = sort(K);  PCH: removed this line.
kmax = K(end);   % PCH: changed "Knew" to "K".

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

% Lower bound(s). Default NaN. Can be specified as either scalar or vector
% of same length as x. If vector, check length/orientation.
if isfield(options,'lbound')
    lbound = options.lbound(:);
    if length(lbound) ~= n && length(lbound) ~= 1
        error('lbound (lower bound) must either be scalar or same length as vector of unknowns x.')
    end
else
    lbound = nan;
end

% Upper bound(s). Default NaN. Can be specified as either scalar or vector
% of same length as x. If vector, check length/orientation.
if isfield(options,'ubound')
    ubound = options.ubound(:);
    if length(ubound) ~= n && length(ubound) ~= 1
        error('ubound (upper bound) must either be scalar or same length as vector of unknowns x.')
    end
else
    ubound = nan;
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

relaxparinput = nan;
if isfield(options,'relaxpar')
    relaxparinput = options.relaxpar;
end

s1 = nan;
if isfield(options,'s1')
     s1 = options.s1;
end

% % If M is given as input
% M = nan;
% if isfield(options,'restart') && isfield(options.restart,'M')
%     M = options.restart.M;
% end
% 
% % If T is given as input
% s = nan;
% if isfield(options,'restart') && isfield(options.restart,'T')
%     s = options.restart.T;
% end

% If weights are given as input
w = nan;
if isfield(options,'w')
    w = options.w;
end


% PCH: I think we must force the user to specify res_dims, otherwise we
% risk that the user forgets to specify it!  (Speking from own experience.)
if strcmp(stoprule,'NCP') && ~isfield(options.stoprule,'res_dims')
    error('options.stoprule.res_dims must be specified for NCP stopping rule.')
end
% res_dims = n;  PCH: removed this line (and it should be "= m").
res_dims = []; % PCH: added this line because "res_dims" must be returned.
if isfield(options,'stoprule') && isfield(options.stoprule,'res_dims')
    res_dims = options.stoprule.res_dims;
    if ~strcmp(stoprule,'NCP')
        % PCH: modified the string below a little bit.
        warning(['options.stoprule.res_dims given but only used by NCP ',...
            'stopping rule, not the one currently specified.']);
    end
end

ncp_smooth = 4;
if isfield(options,'stoprule') && isfield(options.stoprule,'ncp_smooth')
    ncp_smooth = options.stoprule.ncp_smooth;
    if ~strcmp(stoprule,'NCP')
        % PCH: modified the string below a little bit.
        warning(['options.stoprule.ncp_Smooth given but only used by NCP ',...
            'stopping rule, not the one currently specified.']);
    end
end


% Damping.
damp = 0;
if isfield(options,'damping')
    damp = options.damping;
    if damp<0, error('Damping must be positive'), end
end

%% For CART with flagging

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
