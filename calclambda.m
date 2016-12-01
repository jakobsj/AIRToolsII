function [lambda, casel, restart] = calclambda(lambdainput, s1, kmax, atma, n)

% Check if the largest singular value is given.
if isnan(s1)
    % Calculates the largest singular value.
    optionsEIGS.disp = 0;
    sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
else
    sigma1tilde = s1;
end

% Determine the relaxation parameter lambda.
% If lambdainput is nan, set default
if isnan(lambdainput)
    
    % Define a default constant lambda value.
    lambda = 1.9/sigma1tilde^2;
    casel = 1;
    
    % If lambda is a scalar.
elseif ~ischar(lambdainput)
    
    % Checks if the given constant lambde value is unstable.
    if lambdainput <= 0 || lambdainput >= 2/sigma1tilde^2
        warning('MATLAB:UnstableRelaxParam',['The lambda value '...
            'is outside the interval (0,%f)'],2/sigma1tilde^2)
    end
    lambda = lambdainput;
    casel = 1;
else
    % Calculate the lambda value according to the chosen method.
    if strncmpi(lambdainput,'line',4)
        % Method: Line search
        casel = 2;
        
    elseif strncmpi(lambdainput,'psi1',4)
        % Method: ENH psi1.
        casel = 3;
        ite = 0;
        
        % Precalculate the roots.
        z = calczeta(2:kmax-1);
        
        % Define the values for lambda according to the psi1
        % strategy modified or not.
        if strncmpi(lambdainput,'psi1mod',7)
            nu = 2;
            lambda = [sqrt(2); sqrt(2); nu*2*(1-z)]/sigma1tilde^2;
        else
            
            lambda = [sqrt(2); sqrt(2); 2*(1-z)]/sigma1tilde^2;
        end
        
    elseif strncmpi(lambdainput,'psi2',4)
        % Method: ENH psi2.
        casel = 3;
        ite = 0;
        
        % Precalculate the roots.
        kk = 2:kmax-1;
        z = calczeta(kk);
        
        % Define the values for lambda according to the psi2
        % strategy modified or not.
        if strncmpi(lambdainput,'psi2Mod',7)
            nu = 1.5;
            lambda = [sqrt(2); sqrt(2);
                nu*2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
        else
            lambda = [sqrt(2); sqrt(2);
                2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
        end
        
    else
        error(['The chosen relaxation strategy is not valid '...
            'for this method.'])
    end % end check of the class of lambda.
    
end % end check of lambda strategies.

% Save sigma1tilde in restart-struct
restart.s1 = sigma1tilde;