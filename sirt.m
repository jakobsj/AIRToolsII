function [X,info,restart] = sirt(sirt_method, varargin)

% Parse inputs.
% PCH: removed "knew" from output.
% [Afun,b,m,n,K,Knew,kmax,x0,nonneg,boxcon,L,stoprule,taudelta,...
[Afun,b,m,n,K,kmax,x0,nonneg,boxcon,L,stoprule,taudelta,...
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

X = zeros(n,length(K));   % PCH: changed "Knew" to "K".

% TODO Should there be both Knew and K, or reduce to just Knew?
% PCH: fixed this, removed Knew.

% PCH: changed variable names for better consistency:
%      rk -> rkm1 (residual correcponding to iteration k-1)
%      rxk -> rk  (residual corresponding to iteration k)

rk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
[k,K,rkm1,dk] = init_stoprules(stoprule,rk,K,ncp_smooth);

% Do initial check of stopping criteria - probably lambda should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, lambdainput, taudelta, k, kmax, rkm1, dk, res_dims);

% TODO If is aborting here, is output X set? Perhaps make sure x0 is always
% written to the first column of X.

% Calculate lambda and restart
atma = @(x) sfun( Afun( Mfun(Afun(x,'notransp')) , 'transp' ));
[lambda, casel, restart] = calclambda(lambdainput, s1, K, atma, n);

restart.M = M;

% Initialize the values.
xk = x0;
l = 1;   %  PCH: changed "l = 0" to "l = 1".
klast = 0;

while ~stop
    % Update the iteration number k.
    k = k + 1;
    % PCH: deleted the following 3 lines because I redefined the
    % stopping rules ME and NCP, such that we don't need xkm1.
%     if strncmpi(stoprule,'ME',2)  || strncmpi(stoprule,'NC',2)
%         xkm1 = xk;
%     end
    % Compute the current iteration.
    if casel == 1
        % Landweber using constant value of lambda.
        lambdacur = lambda;
        %xk = xk + lambdacur*(Afun(rxk,'transp'));
        % DROP using constant value of lambda.
        %xk = xk + lambda*(s.*(A'*(M.*rxk)));
        xk = xk + lambdacur*(sfun(Afun(Mfun(rk),'transp')));  %(A'*(M.*rxk)));
    elseif casel == 2
        % Landweber using line search.
        %ATrk = Afun(rxk,'transp');
        %lambda = norm(rxk)^2/norm(ATrk)^2;
        %lambdacur = lambda;
        %xk = xk + lambda*ATrk;
        
        % DROP using Line search
        Mrk = Mfun(rk);
        ATMrk = Afun(Mrk,'transp');  %A'*Mrk;        
        ATMrkS = sum(sfun(ATMrk.^2));
        
        lambdacur = (rk'*Mrk)/ATMrkS;
        xk = xk + lambdacur*(sfun(ATMrk));
        
    elseif casel == 3
        % Landweber using psi1 or psi2.
        lambdacur = lambda(k);
        %xk = xk + lambdacur*(Afun(rxk,'transp'));

        %xk = xk + lambdak(ite)*(s.*(A'*(M.*rxk)));
        xk = xk + lambdacur*(sfun(Afun(Mfun(rk),'transp')));  %(A'*(M.*rxk)));
    end % end the different cases of lambda strategies.
    
    % Nonnegativity and box constraints.
    if nonneg, xk = max(xk,0); end
    if boxcon, xk = min(xk,L); end
    
    % New residual.
    rk = b - Afun(xk,'notransp');
    
    % Check stopping rules:
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, lambdacur, taudelta, k, kmax, rkm1, dk, res_dims);
        
    % If the current iteration is requested saved.
    if k == K(l) || stop   % PCH: changed "Knew(l+1)" to "K(l)".
        % l = l + 1;          % PCH: moved this line down.
        % Save the current iteration.
        % PCH rewrote the lines below.
        if strncmpi(stoprule,'ME',2) && stop && info(1)==3
                X(:,l) = xk;   % PCH: "l" ok; changed "xkm1" to "xk".
        elseif strncmpi(stoprule,'NC',2)
            if ~(stop && klast == k-1)
                X(:,l) = xk;   % PCH: "l" ok; changed "xkm1" to "xk".
            else
                % Since the iteration was not saved.
                % PCH: it is a big mystery why this is here!
                l = l - 1;
            end
        else
            X(:,l) = xk;   % PCH: ok
        end
        klast = k;
        l = l + 1;   % PCH: moved line from above to here.
        
    end
end

% Return only the saved iterations.   % PCH: added this comment.
X = X(:,1:l-1);   % PCH: "l-1" because "l" points to next candidate.