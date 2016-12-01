function [X,info,restart] = sirt(sirt_method, varargin)

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,nonneg,boxcon,L,stoprule,taudelta, ...
    lambdainput,s1,M,w,s,res_dims,ncp_smooth] = check_inputs(varargin{:});

% Extract the Mfun and sfun characterizing each SIRT-type method.
if ischar(sirt_method)
    [Mfun,sfun] = get_Mfun_sfun(sirt_method,varargin{1},m,M,w,s);
else
    % Possible to pass in custom SIRT method given by 2-element cell array
    % holding function handles to Mfun and sfun.
    Mfun = sirt_method{1};
    sfun = sirt_method{2};
end

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
[k,K,rkm1,dk] = init_stoprules(stoprule,rk,K,ncp_smooth);

% Do initial check of stopping criteria - probably lambda should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, lambdainput, taudelta, k, kmax, rkm1, dk, res_dims);

% Calculate lambda and restart.
% TODO - NOT NECESSARY FOR SART.
atma = @(x) sfun( Afun( Mfun(Afun(x,'notransp')) , 'transp' ));
[lambda, casel, restart] = calclambda(lambdainput, s1, K, atma, n);

% TODO - how to store M and T/s in third output struct.
restart.M = M;

% Initialize the values.
xk = x0; 
l = 1;   % Pointing to the next iterate number in K to be saved.

while ~stop
    
    % Update the iteration number k.
    k = k + 1;

    % Compute the current iteration depending on lambda strategy.
    if casel == 1
        % SIRT using constant value of lambda.
        lambdacur = lambda;
        xk = xk + lambdacur*(sfun(Afun(Mfun(rk),'transp')));
        
    elseif casel == 2
        % SIRT using line search.
        Mrk = Mfun(rk);
        ATMrk = Afun(Mrk,'transp');  %A'*Mrk;        
        ATMrkS = sum(sfun(ATMrk.^2));
        lambdacur = (rk'*Mrk)/ATMrkS;
        xk = xk + lambdacur*(sfun(ATMrk));
        
    elseif casel == 3
        % SIRT using psi1 or psi2.
        lambdacur = lambda(k);
        xk = xk + lambdacur*(sfun(Afun(Mfun(rk),'transp')));
    end % end the different cases of lambda strategies.
    
    % Nonnegativity and box constraints.
    if nonneg, xk = max(xk,0); end
    if boxcon, xk = min(xk,L); end
    
    % New residual.
    rk = b - Afun(xk,'notransp');
    
    % Check stopping rules.
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, lambdacur, taudelta, k, kmax, rkm1, dk, res_dims);
        
    % If the current iteration is requested saved.
    if k == K(l) || stop
        X(:,l) = xk;
        l = l + 1;
    end
end

% Return only the saved iterations: Only to "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);