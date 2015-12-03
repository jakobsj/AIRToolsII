function tau = trainDPME(A,b_exact,x_exact,method,type,delta,s,options)
%TRAINDPME Training method for the stopping rules DP and ME
%
%   tau = trainDPME(A,b_exact,x_exact,method,type,delta,s)
%   tau = trainDPME(A,b_exact,x_exact,method,type,delta,s,options)
%
% This function determines the parameter tau for a given method, for
% the stopping rule DP or ME, by training on a problem given by A, b_exact,
% and x_exact. The noise level is delta, s samples of the noisy right-hand
% side are created, and from these samples the value of tau is determined.
%
% Input:
%   A           m times n matrix.
%   b_exact     m times 1 vector containing the exact rhs.
%   x_exact     n times 1 vectir containing the exact solution.
%   method      Function handle to a SIRT or ART method.
%   type        String that should be either 'DP' for Discrepancy Principle
%               or 'ME' for the Monotone Error rule.  ME can only be chosen
%               if method is a SIRT-method.
%   delta       Scalar denoting the noise level || e ||_2.
%   s           The number of used samples.
%   options     Struct used in the call of the method.  For this strategy
%               the field stoprule cannot be used.
%
% Output:
%   tau         Scalar containing the trained parameter.
%
% See also: cav, cimmino, drop, kaczmarz, landweber, randkaczmarz, sart,
% symkaczmarz.

% Maria Saxild-Hansen and Per Chr. Hansen, June 23, 2010, DTU Compute.
% Bug fix by Camilla Trinderup.

% Reference: T. Elfving and T. Nikazad, Stopping rules for Landweber-type
% iteration, Inverse Problems, 23 (2007), pp. 1417-1432.

if nargin == 8
    if isfield(options,'stoprule')
        options = rmfield(options,'stoprule');
    end
end

m = size(A,1);

% Generate s noisy samples of the rhs.
% Debugging
e = randn(m,s);
e = delta*e./(ones(m,1)*sqrt(sum(e.*e)));
b = b_exact*ones(1,s) + e;

% The norm of the exact solution.
normx_ex = norm(x_exact);

% Initialize the tauh vector.
tauh = zeros(s,1);

% Handle the different methods according to their type.
switch func2str(method)
    case {'landweber','cimmino','cav','drop','sart'}
        
        % Determine the M12 and its norm.
        [ignore ignore restart] = method(A,b(:,1),1);
        M = restart.M;
        if ~isempty(M)
            M12 = sqrt(M);
            M12norm = max(M12);
        else
            M12 = 1;
            M12norm = 1;
        end
        options.restart = restart;
        options.stoprule.type = 'none';
        % The maximum number of iterations for SIRT methods.
        kmax = 1000;
        
    case {'kaczmarz','symkaczmarz','randkaczmarz'}
        
        % Check that the ART methods only use DP.
        if ~strncmpi(type,'DP',2)
            error('For this method you can only use the type DP')
        else
            M12 = 1;
            M12norm = 1;
        end
        options.stoprule.type = 'none';
        % The maximum number of iterations for ART methods.
        kmax = 100;
        
    otherwise
        error('Unknown method.')
        
end

% The maximum number of iterations.
if kmax < 10
    kstep = kmax;
else
    kstep = 10;
end

% Loop over all generated samples of the rhs b.
for i = 1:s
    
    % Prepare for solving with the i'th rhs.
    bi = b(:,i);
    X = zeros(size(x_exact,1),1);
    x0 = zeros(size(x_exact));
    k = kstep;
    K = 1:k;
    Ek = [];
    
    stop = 0;
    
    % As long as the minimum residual is not found.
    while ~stop
        
        % Perform the K iterations.
        Xny = method(A,bi,K,x0,options);
        
        % Calculate the relative error.
        xd = Xny - repmat(x_exact,1,size(Xny,2));
        Ek = [Ek sqrt(sum(xd.*xd,1))/normx_ex];
        X = [X Xny];
        
        % Find the minimum relative error.
        [minEk I] = min(Ek);
        
        % If the minimum error is not in the last iteration, then the
        % loop is stopped.
        if I ~= k || k >= kmax
            
            % Stop the loop and define the optimal number of iterations.
            stop = 1;
            kopt = I;
            
        else
            
            % Continue the iterations, and update iterations and starting
            % value.
            if k+kstep > kmax
                K = 1:(kmax-k);
                k = kmax;
            else
                k = k + kstep;
                % CT's addition
                
            end
            x0 = Xny(:,end);
            
        end
    end
    
    % Calculate the residual for the optimal iteration and the previous one.
    rk = M12.*(bi-A*X(:,kopt+1));
    rkm1 = M12.*(bi-A*X(:,kopt));
    
    % Split in the types DP and ME.
    if strncmpi(type,'DP',2)
        % Real formula is d(k)/(delta*norm(M^{1/2})*norm(rk)), but since
        % dk = norm(rk)^2 for DP we get:
        
        Rk = norm(rk)/(delta*M12norm);
        Rkm1 = norm(rkm1)/(delta*M12norm);
        
    elseif strncmpi(type,'ME',2)
        
        % Ensuring that there is enough iterations to find the next
        % residual
        if kopt +2 > size(X,2)
            Xny = method(A,bi,2,Xny(:,end),options);
            X = [X Xny]; 
        end
        
        % Define the residual for the iteration after the optimal.
        rkp1 = M12.*(bi-A*X(:,kopt+2));
               
        Rk = (rk'*(rk+rkp1))/(delta*M12norm*norm(rk));
        Rkm1 = (rkm1'*(rkm1+rk))/(delta*M12norm*norm(rkm1));
        
    else
        error('The chosen type is invalid')
    end
    
    
    % Define the tau value for the current right-hand side.
    tauh(i) = mean([Rk Rkm1]);
end

% Determine the optimal tau value.
tau = mean(tauh);