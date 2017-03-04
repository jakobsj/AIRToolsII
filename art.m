function [X,info] = art(art_method, varargin)

% Set default ART method to be kaczmarz.
if isempty(art_method)
    art_method = 'kaczmarz';
end

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    s1,w,res_dims,rkm1,dk,damp] = check_inputs(varargin{:});

% Special check for symkaczmarz: number of iterations must be even
if ischar(art_method) && strncmpi(art_method,'sym',3)
    if any(mod(K,2))
        error('For symkaczmarz only even iteration numbers can be requested.');
    else
        % Since in symkaczmarz a sweep is top->bottom->top it is twice as
        % expensive as other methods, and one sweep is counted as two
        % iterations. For code reasons, we do this by running half the
        % number of iterations, then later doubling iteration numbers
        % again.
        K = K/2;
        kmax = kmax/2;
    end    
end

% If A is matrix, it is more efficient for ART to work with transposed A.
A = varargin{1};
if ~isa(A,'function_handle')
    A = A';
end

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialization before iterations.
k = 0;   % Iteration counter
xk = x0; % Use initial vector.
l = 1;   % Pointing to the next iterate number in K to be saved.

% Do initial check of stopping criteria - probably relaxpar should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Compute the relaxation parameter to be used throughout iterations.
relaxpar = calc_relaxpar(relaxparinput);

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
is_randkaczmarz = false;
if ischar(art_method)
    switch lower(art_method)
        case 'kaczmarz'
            I = find(normAi>0);
        case 'symkaczmarz'
            I = find(normAi>0);
            I = [I, I(end-1:-1:2)];
        case 'randkaczmarz'
            is_randkaczmarz = true;
            I = find(normAi>0);
        otherwise
            error('Unknown ART method specified')
    end
else
    % Custom ART method specified by user giving the row order I as
    % first input instead of the string with a particular ART method name.
    I = art_method(:)';
    I = I(normAi(I)>0);
end

% Apply damping.
normAi = normAi + damp*max(normAi);

% Special for randkaczmarz - set up random row selection 
if is_randkaczmarz
    cumul = cumsum(normAi/sum(normAi));
    
    % If all rows have approximately the same norm, treat them as having the
    % same norm when computing the random row index (which is much faster).
    if norm(cumul-(1:m)/m,inf) < 0.05  % Arbitrary threshold.
        fast = true;
    else
        fast = false;
    end
end    

% Main ART loop.
while ~stop
    
    % Update the iteration number k.
    k = k + 1;
    
    % Special for randkaczmarz - need random permutation or rows.
    if is_randkaczmarz
        I_torun = I(randperm(length(I)));
    else
        I_torun = I;
    end
    
    % The Kaczmarz sweep.
    for i = I_torun
        
        % Special for randkaczmarz - 
        if is_randkaczmarz
            % The random row index.
            if fast
                ri = i;
            else
                ri = sum(cumul<rand)+1;
            end
        else
            % All other methods, just use i.
            ri = i;
        end
        
        if isa(A,'function_handle')
            e = zeros(m,1); e(ri) = 1;
            ai = Afun(e,'transp');  % ai is a column vector.
        else
            ai = A(:,ri); % Remember that A is transposed.
        end
        xk = xk + (relaxpar*(b(ri) - ai'*xk)/normAi(ri))*ai;
        
        % Enforce any lower and upper bounds (scalars or xk-sized vectors)
        if ~isempty(lbound)
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

% Special for symkaczmarz: double finaliter
if ischar(art_method) && strncmpi(art_method,'sym',3)
    info.finaliter = info.finaliter*2;
end

% List of iterates saved: all in K smaller than the final, and the final.
info.itersaved = [K(K<info.finaliter), info.finaliter];
