function [relaxpar, casel, sigma1tilde] = calcrelaxpar_sirt(relaxparinput, s1, kmax, atma, n)

% Check if the largest singular value is given.
if isnan(s1)
    % Calculates the largest singular value.
    optionsEIGS.disp = 0;
    sigma1tilde = sqrt( eigs(atma,n,1,'lm',optionsEIGS) );
else
    sigma1tilde = s1;
end

% Determine the relaxation parameter relaxpar.
% If relaxparinput is nan, set default
if isnan(relaxparinput)
    
    % Define a default constant relaxpar value.
    relaxpar = 1.9/sigma1tilde^2;
    casel = 1;
    
    % If relaxpar is a scalar.
elseif ~ischar(relaxparinput)
    
    % Checks if the given constant lambde value is unstable.
    if relaxparinput <= 0 || relaxparinput >= 2/sigma1tilde^2
        warning('MATLAB:UnstableRelaxParam',['The relaxpar value '...
            'is outside the interval (0,%f)'],2/sigma1tilde^2)
    end
    relaxpar = relaxparinput;
    casel = 1;
else
    % Calculate the relaxpar value according to the chosen method.
    if strncmpi(relaxparinput,'line',4)
        % Method: Line search
        casel = 2;
        
    elseif strncmpi(relaxparinput,'psi1',4)
        % Method: ENH psi1.
        casel = 3;
        ite = 0;
        
        % Precalculate the roots.
        z = calczeta(2:kmax-1);
        
        % Define the values for relaxpar according to the psi1
        % strategy modified or not.
        if strncmpi(relaxparinput,'psi1mod',7)
            nu = 2;
            relaxpar = [sqrt(2); sqrt(2); nu*2*(1-z)]/sigma1tilde^2;
        else
            
            relaxpar = [sqrt(2); sqrt(2); 2*(1-z)]/sigma1tilde^2;
        end
        
    elseif strncmpi(relaxparinput,'psi2',4)
        % Method: ENH psi2.
        casel = 3;
        ite = 0;
        
        % Precalculate the roots.
        kk = 2:kmax-1;
        z = calczeta(kk);
        
        % Define the values for relaxpar according to the psi2
        % strategy modified or not.
        if strncmpi(relaxparinput,'psi2Mod',7)
            nu = 1.5;
            relaxpar = [sqrt(2); sqrt(2);
                nu*2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
        else
            relaxpar = [sqrt(2); sqrt(2);
                2*(1-z)./((1-z.^(kk')).^2)]/sigma1tilde^2;
        end
        
    else
        error(['The chosen relaxation strategy is not valid '...
            'for this method.'])
    end % end check of the class of relaxpar.
    
end % end check of relaxpar strategies.
