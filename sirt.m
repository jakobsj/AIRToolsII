function [X,info,restart] = sirt(sirt_method, varargin)




% Parse inputs.
[Afun,b,m,n,K,Knew,kmax,x0,nonneg,boxcon,L,stoprule,taudelta,...
    lambdainput,s1,M,w,s,res_dims,ncp_smooth] = check_inputs(varargin{:});

% Extract the Mfun and sfun characterizing each SIRT-type method
if ischar(sirt_method)
    [Mfun,sfun] = get_Mfun_sfun(sirt_method,varargin{1},m,M,w,s);
else
    % Possible to pass in custom SIRT method given by 2-element cell array
    % holding function handles to Mfun and sfun.
    Mfun = sirt_method{1};
    sfun = sirt_method{2};
end


X = zeros(n,length(Knew));

% TODO Should there be both Knew and K, or reduce to just Knew?

rxk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
[k,K,rk,dk] = init_stoprules(stoprule,rxk,K,ncp_smooth);

% Do initial check of stopping criteria - probably lambda should be set
% before this, perhaps just to nan.
[stop, info, rk, dk] = check_stoprules(...
    stoprule, rxk, lambdainput, taudelta, k, kmax, rk, dk, res_dims);

% TODO If is aborting here, is output X set? Perhaps make sure x0 is always
% written to the first column of X.

% Calculate lambda and restart
atma = @(x) sfun( Afun( Mfun(Afun(x,'notransp')) , 'transp' ));
[lambda, casel, restart] = calclambda(lambdainput, s1, K, atma, n);

restart.M = M;

% Initialize the values.
xk = x0;
l = 0;
klast = 0;

while ~stop
    % Update the iteration number k.
    k = k + 1;
    if strncmpi(stoprule,'ME',2)  || strncmpi(stoprule,'NC',2)
        xkm1 = xk;
    end
    % Compute the current iteration.
    if casel == 1
        % Landweber using constant value of lambda.
        lambdacur = lambda;
        %xk = xk + lambdacur*(Afun(rxk,'transp'));
        % DROP using constant value of lambda.
        %xk = xk + lambda*(s.*(A'*(M.*rxk)));
        xk = xk + lambdacur*(sfun(Afun(Mfun(rxk),'transp')));  %(A'*(M.*rxk)));
    elseif casel == 2
        % Landweber using line search.
        %ATrk = Afun(rxk,'transp');
        %lambda = norm(rxk)^2/norm(ATrk)^2;
        %lambdacur = lambda;
        %xk = xk + lambda*ATrk;
        
        % DROP using Line search
        Mrk = Mfun(rxk);
        ATMrk = Afun(Mrk,'transp');  %A'*Mrk;        
        ATMrkS = sum(sfun(ATMrk.^2));
        
        lambdacur = (rxk'*Mrk)/ATMrkS;
        xk = xk + lambdacur*(sfun(ATMrk));
        
    elseif casel == 3
        % Landweber using psi1 or psi2.
        lambdacur = lambda(k);
        %xk = xk + lambdacur*(Afun(rxk,'transp'));

        %xk = xk + lambdak(ite)*(s.*(A'*(M.*rxk)));
        xk = xk + lambdacur*(sfun(Afun(Mfun(rxk),'transp')));  %(A'*(M.*rxk)));
    end % end the different cases of lambda strategies.
    
    % Nonnegativity and box constraints.
    if nonneg, xk = max(xk,0); end
    if boxcon, xk = min(xk,L); end
    
    % New residual.
    rxk = b - Afun(xk,'notransp');
    
    % Check stopping rules:
    [stop,info,rk,dk] = check_stoprules(...
        stoprule, rxk, lambdacur, taudelta, k, kmax, rk, dk, res_dims);
        
    % If the current iteration is requested saved.
    if k == Knew(l+1) || stop
        l = l + 1;
        % Save the current iteration.
        % PCH rewrote the lines below.
        if strncmpi(stoprule,'ME',2) && stop && info(1)==3
                X(:,l) = xkm1;
        elseif strncmpi(stoprule,'NC',2)
            if ~(stop && klast == k-1)
                X(:,l) = xkm1;
            else
                % Since the iteration was not saved.
                l = l - 1;
            end
        else
            X(:,l) = xk;
        end
        klast = k;
        
    end
end
X = X(:,1:l);