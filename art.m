function [X,info] = art(art_method, varargin)

% Set default ART method to be kaczmarz.
if isempty(art_method)
    art_method = 'kaczmarz';
end

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    s1,w,res_dims,ncp_smooth,damp] = check_inputs(varargin{:});

A = varargin{1};
if ~isa(A,'function_handle')
    A = A';
end

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
if strcmpi(stoprule,'ME')
    error('Stopping rule ME is not available for ART methods.')
end
[k,K,rkm1,dk] = init_stoprules(stoprule,rk,K,ncp_smooth);

% Do initial check of stopping criteria - probably relaxpar should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Compute the relaxation parameter to be used throughout iterations.
relaxpar = calcrelaxpar_art(relaxparinput);

% Calculate the norm of each row in A. This calculation can require a
% lot of memory.
if ~isa(A,'function_handle')
    normAi = full(abs(sum(A.*A,1)));  % Remember that A is transposed.
else
    normAi = zeros(1,m);
    for i = 1:m
        e = zeros(m,1);
        e(i) = 1;
        v = A(e,'transp');
        normAi(i) = norm(v)^2;
    end
end

% Depending on ART method, set the row order.
if ischar(art_method)
    switch lower(art_method)
        case 'kaczmarz'
            I = find(normAi>0);
        case 'symkaczmarz'
            I = find(normAi>0);
            I = [I, I(end-1:-1:2)];
        case 'randkaczmarz'
            error('not implemented')
        otherwise
            error('Unknown ART method specified')
    end
else
    % Custom ART method specified by user giving the row order I as
    % first input instead of the string with a particular ART method name.
    I = art_method;
end

% Apply damping.
normAi = normAi + damp*max(normAi);

% Initialization before iterations.
xk = x0;
l = 1;

while ~stop
    
    % Update the iteration number k.
    k = k + 1;
    
    % The Kaczmarz sweep.
    for i = I
        if isa(A,'function_handle')
            e = zeros(m,1); e(i) = 1;
            ai = Afun(e,'transp');  % ai is a column vector.
        else
            ai = A(:,i); % Remember that A is transposed.
        end
        xk = xk + (relaxpar*(b(i) - ai'*xk)/normAi(i))*ai;
        
        % Enforce any lower and upper bounds (scalars or xk-sized vectors)
        if ~isnan(lbound)
            xk = max(xk,lbound);
        end
        if ~isnan(ubound)
            xk = min(xk,ubound);
        end
    end
    
    % New residual.
    rk = b - Afun(xk,'notransp');
    
    % Check stopping rules.
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, relaxpar, taudelta, k, kmax, rkm1, dk, res_dims);
    
    % If the current iteration is requested saved.
    if k == K(l) || stop
        X(:,l) = xk;
        l = l + 1;
    end
end

% Return only the saved iterations: Only to "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);
