function [X,info,ext_info] = sirt(sirt_method, varargin)

% options.s1 now instead of options.restart.s1: largest singval

% Set default SIRT method to be sart.
if isempty(sirt_method)
    sirt_method = 'sart';
end

% Parse inputs.
[Afun,b,m,n,K,kmax,x0,lbound,ubound,stoprule,taudelta, relaxparinput, ...
    s1,w,res_dims,ncp_smooth] = check_inputs(varargin{:});

% Extract the Mfun and sfun characterizing each SIRT-type method.
if ischar(sirt_method)
    [Mfun,Tfun] = get_mfun_tfun(sirt_method,varargin{1},m,n,w);
else
    % Possible to pass in custom SIRT method given by 2-element cell array
    % holding function handles to Mfun and sfun instead of string input.
    Mfun = sirt_method{1};
    Tfun = sirt_method{2};
end

% Initialize array to hold requested iterates.
X = zeros(n,length(K));

% Residual of initial guess.
rk = b - Afun(x0,'notransp');

% Initialize for stopping rules.
[k,rkm1,dk] = init_stoprules(stoprule,rk,ncp_smooth);

% Do initial check of stopping criteria - probably relaxpar should be set
% before this, perhaps just to nan.
[stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxparinput, taudelta, k, kmax, rkm1, dk, res_dims);

% Calculate relaxpar. Special case for SART largest singval is 1.
atma = @(x) Tfun( Afun( Mfun(Afun(x,'notransp')) , 'transp' ));
if strcmpi(sirt_method,'sart')
    s1 = 1;
end
[relaxpar, casel, sigma1tilde] = calc_relaxpar_sirt(relaxparinput, s1, kmax, atma, n);

% Store M and T in third output struct.
if nargout > 2
    ext_info.M = Mfun(ones(m,1));
    ext_info.T = Tfun(ones(n,1));
end

% Initialize the values.
xk = x0; 
l = 1;   % Pointing to the next iterate number in K to be saved.

while ~stop
    
    % Update the iteration number k.
    k = k + 1;

    % Compute the current iteration depending on relaxpar strategy.
    Mrk = Mfun(rk);
    ATMrk = Afun(Mrk,'transp');  %A'*Mrk;  
    if casel == 1
        % SIRT using constant value of relaxpar.
        relaxparcur = relaxpar;        
    elseif casel == 2
        % SIRT using line search.    
        ATMrkS = sum(Tfun(ATMrk.^2));
        relaxparcur = (rk'*Mrk)/ATMrkS;
    elseif casel == 3
        % SIRT using psi1 or psi2.
        relaxparcur = relaxpar(k);
    end % end the different cases of relaxpar strategies.
    
    % The update step with current relaxpar
    xk = xk + relaxparcur*(Tfun(ATMrk));
    
    % Enforce any lower and upper bounds (scalars or xk-sized vectors)
    if ~isempty(lbound)
        xk = max(xk,lbound);
    end
    if ~isempty(ubound)
        xk = min(xk,ubound);
    end
    
    % New residual.
    rk = b - Afun(xk,'notransp');
    
    % Check stopping rules.
    [stop,info,rkm1,dk] = check_stoprules(...
        stoprule, rk, relaxparcur, taudelta, k, kmax, rkm1, dk, res_dims);
        
    % If the current iteration is requested saved.
    if k == K(l) || stop
        X(:,l) = xk;
        l = l + 1;
    end
end

% Return only the saved iterations: Only to "l-1" because "l" now points to
% next candidate.
X = X(:,1:l-1);


% Save to info:

% Largest singular value determined
info.s1 = sigma1tilde;

% List of iterates saved: all in K smaller than the final, and the final.
info.itersaved = [K(K<info.finaliter), info.finaliter];
