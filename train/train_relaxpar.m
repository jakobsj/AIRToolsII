function relaxpar = train_relaxpar(A,b,x_ex,method,kmax,options)
%TRAIN_RELAXPAR  Training to find optimal relaxpar for ART/CART/SIRT methods
%
%   relaxpar = train_relaxpar(A,b,x_ex,method,kmax)
%   relaxpar = train_relaxpar(A,b,x_ex,method,kmax,options)
%
% This function determines the optimal value of the relaxation parameter
% relaxpar for one of the SIRT methods cav, cimmino, drop, landweber, and
% sart, one of the ART methods kaczmarz, symkaczmarz and randkaczmarz, and
% the CART  method columnaction. The optimal value of relaxpar is defined
% as the value that gives rise to the fastest convergence to the smallest
% error  in the solution.
%
% Input:
%   A           m times n matrix or function handle to matrix-free version.
%   b           m times 1 vector containing the right-hand side.
%   x_ex        n times 1 vector containing the exact solution.
%   method      Function handle to one of the SIRT/ART/CART methods.
%   kmax        Scalar that determines the maximum number of iterations
%               of the used SIRT/ART/CART-method.
%   options     Struct used in the call of the method. For this strategy
%               the fields stoprule and relaxpar cannot be used.
%
% Output:   
%   relaxpar    Scalar containing the found relaxpar value.
%
% See also: demo_training, train_dpme.

% Reference: P. C. Hansen and M. Saxild-Hansen, AIR Tools - A MATLAB package
% of algebraic iterative reconstruction methods, J. Comp. Appl. Math., 236
% (2012), pp. 2167-2178.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% Input check: ensure that no stoprule is chosen.
if nargin == 6
    options.stoprule.type = 'none';
    if isfield(options,'relaxpar')
        options = rmfield(options,'relaxpar');
    end
else
    options = [];
end 

% Fudge factor.
fudge = 0.02;

% Prepare for step 1.
stepk = 10;
K = 1:stepk;
k = 0;
stop = 0;
Er = [];
x0 = [];

% Step 1: Determine the minimum error for the default relaxpar.
while ~stop
    % Ensure only kmax iterations is performed.
    if k+stepk > kmax
        K = 1:(kmax-k);
        k = kmax;
    else
        % Update the number of iterations.
        k = k + stepk;
    end
    
    % If the first run in the loop.
    if k <= stepk
        switch func2str(method)
            case {'landweber','cimmino','cav','drop','sart'}
                
                % Compute the solutions to the first K values.
                [Xnew, info] = method(A,b,K,x0,options);
                
                % Assign the computed rho to options, such that rho does
                % not require recalculation in the next call of the method.
                options.rho = info.rho;
                relaxparmax = 2/info.rho;
                
                is_sirt = true;
            
            case {'columnaction','kaczmarz','randkaczmarz','symkaczmarz'}
                
                % Compute the solutions for the first iterations.
                options.relaxpar = 0.25;
                Xnew = method(A,b,K,x0,options);
                relaxparmax = 2;
                
                is_sirt = false;
                
            otherwise
                error(['Unknown method ',func2str(method)])
        end
    else
        % Compute the solutions to the next K values.
        Xnew = method(A,b,K,x0,options);
    end
    
    % Compute the errors for the new solutions and collect all the errors.
    deltax = Xnew - repmat(x_ex,1,size(Xnew,2));
    Er = [Er, sqrt(sum(deltax.*deltax,1))];
    
    % Find the minimum error.
    [minE, kopt] = min(Er);
    
    % If the minimum error is in the last iteration, then repeat the loop
    % with stepk additional iterations.  Otherwise the minimum error is
    % found and the loop is terminated.
    if kopt == k && k < kmax
        x0 = Xnew(:,end);
    else
        stop = 1;
    end
end

% Minimum error upper bound.
meub = (1+fudge)*minE;

% Step 2: Find the relaxpar for which the minimum error is below meub
% and kopt is minimum.

% Define the search interval for relaxpar.
r = (3-sqrt(5))/2;
alpha = 0;
beta = relaxparmax;
alphap = alpha + r*(beta-alpha);
betap = alpha + (1-r)*(beta-alpha);

% Calculate err_alphap and k_alphap.
options.relaxpar = alphap;
xnew = method(A,b,1:kopt,[],options);
ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
k_alphap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
err_alphap = ee(k_alphap);

% Calculate err_betap and k_betap.
options.relaxpar = betap;
xnew = method(A,b,1:kopt,[],options);
ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
k_betap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
err_betap = ee(k_betap);

% Augmented version of golden section search.  For the ART and CART methods
% we favor smaller relaxation parameters because they tend to give faster
% semi-convergence; hence this order of the first two if-elseif blocks. 
% For the SIRT methods we favor larger relaxation parameters, causing the
% opposite order of the first two if-elseif blocks, while the last two have
% same order.
while abs(alphap-betap) > relaxparmax*0.01
    
    if is_sirt
        if err_alphap > meub
            
            % [alpha alphap betap beta] <---- [alphap betap z beta]
            z = alphap + (1-r)*(abs(beta-alphap));
            alpha = alphap;  alphap = betap;  betap = z;
            k_alphap = k_betap;  err_alphap = err_betap;
            
            options.relaxpar = betap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_betap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_betap = ee(k_betap);
            
        elseif err_betap > meub
            
            % [alpha alphap betap beta] <---- [alpha z alphap betap]
            z = alpha + r*(abs(betap-alpha));
            beta = betap;  betap = alphap;  alphap = z;
            k_betap = k_alphap;  err_betap = err_alphap;
            
            options.relaxpar = alphap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_alphap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_alphap = ee(k_alphap);
            
        elseif k_alphap >= k_betap
            
            % [alpha alphap betap beta] <---- [alphap betap z beta]
            z = alphap + (1-r)*(abs(beta-alphap));
            alpha = alphap;  alphap = betap;  betap = z;
            k_alphap = k_betap;  err_alphap = err_betap;
            
            options.relaxpar = betap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_betap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_betap = ee(k_betap);
            
        else
            
            % [alpha alphap betap beta] <---- [alpha z alphap betap]
            z = alpha + r*(abs(betap-alpha));
            beta = betap;  betap = alphap;  alphap = z;
            k_betap = k_alphap;  err_betap = err_alphap;
            
            options.relaxpar = alphap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_alphap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_alphap = ee(k_alphap);
            
        end
    else  % Then it is ART/CART
        if err_betap > meub
            
            % [alpha alphap betap beta] <---- [alpha z alphap betap]
            z = alpha + r*(abs(betap-alpha));
            beta = betap;  betap = alphap;  alphap = z;
            k_betap = k_alphap;  err_betap = err_alphap;
            
            options.relaxpar = alphap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_alphap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_alphap = ee(k_alphap);
            
        elseif err_alphap > meub
            
            % [alpha alphap betap beta] <---- [alphap betap z beta]
            z = alphap + (1-r)*(abs(beta-alphap));
            alpha = alphap;  alphap = betap;  betap = z;
            k_alphap = k_betap;  err_alphap = err_betap;
            
            options.relaxpar = betap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_betap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_betap = ee(k_betap);
            
        elseif k_alphap >= k_betap
            
            % [alpha alphap betap beta] <---- [alphap betap z beta]
            z = alphap + (1-r)*(abs(beta-alphap));
            alpha = alphap;  alphap = betap;  betap = z;
            k_alphap = k_betap;  err_alphap = err_betap;
            
            options.relaxpar = betap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_betap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_betap = ee(k_betap);
            
        else
            
            % [alpha alphap betap beta] <---- [alpha z alphap betap]
            z = alpha + r*(abs(betap-alpha));
            beta = betap;  betap = alphap;  alphap = z;
            k_betap = k_alphap;  err_betap = err_alphap;
            
            options.relaxpar = alphap;
            xnew = method(A,b,1:kopt,[],options);
            ee = sqrt(sum((xnew-repmat(x_ex,1,kopt)).^2,1));
            k_alphap = sum(ee(1:end-1)>meub & diff(ee)<0) + 1;
            err_alphap = ee(k_alphap);
            
        end
    end
end

% Define the optimal relaxpar.
relaxpar = (alphap+betap)/2;