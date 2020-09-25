function [relaxpar, casel, rho] = ...
    calc_relaxpar(relaxparinput, rho, kmax, atma, n)
%CALC_RELAXPAR  Compute the relaxation parameter for all methods
%
%   relaxpar = calc_relaxpar(relaxparinput)
%   [relaxpar, casel, rho] = ...
%              calc_relaxpar(relaxparinput, rhoinput, kmax, atma, n)
%
% The short form is used by art and cart and sets the relaxpar either to
% default values of 1 or 0.25, respectively, or assigns a value given by
% the user. If the value given by the user is outside the allowed interval,
% a warning is given.
% 
% The long form is used by sirt, and inputs and outputs are explained below.
% 
% Input:
%    relaxparinput  Any relaxation parameter or flag determining method 
%                   to use for determining relaxpar as specified by user.
%    rho            The spectral radius of the iteration matrix, if given
%                   by user.
%    kmax           The maximum number of SIRT iterations.
%    atma           A function handle to the iteration matrix that
%                   characterizes the SIRT method, for which to compute
%                   the spectral radius.
%    n              The number of columns in A.
%    
% Output:
%    relaxpar       The computed relaxation parameter, or a vector of
%                   iteration-dependent relaxation parameters.
%    casel          A flag indicating whether a constant lambda is
%                   returned (casel=1), line search is to be used (casel=2)
%                   or the psi1/psi2 strategies to be used (casel=3).
%    rho            Computed spectral radius of the iteration matrix.
%
% See also: art, cart, sirt

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

% First determine whether called from art, cart or sirt. Just testing for
% the specific method name such as cimmino would not cover the case of
% custom methods.
stack = dbstack;

switch stack(2).name
    
    case 'art'
        % Default choice 1. If user gave as input, use that value and 
        % throw warning if outside [0,2], but proceed.
        if isempty(relaxparinput)
            relaxpar = 1;
        elseif isnumeric(relaxparinput)
            if (relaxparinput <= 0 || relaxparinput >= 2)
                warning('MATLAB:UnstableRelaxParam',...
                    'The relaxpar value is outside the interval (0,2)');
            end
            relaxpar = relaxparinput;
        elseif isa(relaxparinput,'function_handle')
            relaxpar = -1000;  % Place holder; function handle will be used
        else
            error('MATLAB:IllegalRelaxParam',...
                ['For ART methods the relaxpar must be a scalar or ',...
                'function handle.'])
        end
        
    case 'cart'
        % Default choice 0.25. If user gave as input, use that value and
        % throw warning if outside [0,2], but proceed.
        if isempty(relaxparinput)
            relaxpar = 0.25;
        elseif isnumeric(relaxparinput)
            if relaxparinput <= 0 || relaxparinput >= 2
                warning('MATLAB:UnstableRelaxParam',...
                    'The relaxpar value is outside the interval (0,2)');
            end
            relaxpar = relaxparinput;
        else
            error('MATLAB:IllegalRelaxParam',...
                'For CART methods the relaxpar must be a scalar.')
        end
        
    case 'sirt'
        % Check if the spectral radius is given.
        if isnan(rho)
            % If not, calculate the largest singular value.
            optionsEIGS.disp = 0;
            optionsEIGS.tol = 1e-4;
            
            % Save existing random number generator settings to be restored
            % after having specified a fixed random seed.  This ensures a
            % deterministic result (eigs without this is non-deterministic
            % due to starting vector chosen using rand).
            scurr = rng(0);
            rho = eigs(atma,n,1,'lm',optionsEIGS);
            rng(scurr);
        end
        
        % Determine the relaxation parameter relaxpar.
        % If relaxparinput is nan, set default.
        if isempty(relaxparinput)
            
            % Define a default constant relaxpar value.
            relaxpar = 1.9/rho;
            casel = 1;
            
            % If relaxpar is a scalar.
        elseif ~ischar(relaxparinput)
            
            % Checks if the given constant lambde value is illegal.
            if relaxparinput <= 0 || relaxparinput >= 2/rho
                warning('MATLAB:UnstableRelaxParam',...
                    ['The relaxpar value '...
                    'is outside the interval (0,%f)'],2/rho)
            end
            relaxpar = relaxparinput;
            casel = 1;
            
        else
            % Calculate the relaxpar value according to the chosen method.
            if strncmpi(relaxparinput,'line',4)
                % Method: Line search.
                casel = 2;
                relaxpar = [];
                
            elseif strncmpi(relaxparinput,'psi1',4)
                % Method: ENH psi1.
                casel = 3;
                
                % Precalculate the roots.
                z = calczeta(2:kmax-1);
                
                % Define the values for relaxpar according to the psi1
                % strategy modified or not.
                if strncmpi(relaxparinput,'psi1mod',7)
                    nu = 2;
                    relaxpar = [sqrt(2); sqrt(2); nu*2*(1-z)]/rho;
                else
                    relaxpar = [sqrt(2); sqrt(2); 2*(1-z)]/rho;
                end
                
            elseif strncmpi(relaxparinput,'psi2',4)
                % Method: ENH psi2.
                casel = 3;
                
                % Precalculate the roots.
                kk = 2:kmax-1;
                z = calczeta(kk);
                
                % Define the values for relaxpar according to the psi2
                % strategy modified or not.
                if strncmpi(relaxparinput,'psi2Mod',7)
                    nu = 1.5;
                    relaxpar = [sqrt(2); sqrt(2);
                        nu*2*(1-z)./((1-z.^(kk')).^2)]/rho;
                else
                    relaxpar = [sqrt(2); sqrt(2);
                        2*(1-z)./((1-z.^(kk')).^2)]/rho;
                end
                
            else
                error(['The chosen relaxation strategy is not valid '...
                    'for this method.'])
            end % end check of the class of relaxpar.
            
        end % end check of relaxpar strategies.
        
    otherwise
        error(['This aux. function is only to be called ',...
            'from within art, cart or sirt.'])
end

function z = calczeta(k)
%CALCZETA Calculates a specific root of a certain polynomial
%
%   z = calczeta(k)
%
% Calculates the unique root in the interval [0 1] of the polynomial equation
%       (2k -1)z^{k-1} - ( z^{k-2} + ... + z + 1 ) = 0
% by means of Newton's method.
%
% Input:
%   k  either a sorted vector or a scalar that determines the degree 
%      of the polynomial g(z).
% Output:
%   z  the found root(s).
%  
% See also: cav, cimmino, drop, landweber, sart, symkaczmarz.

% Maria Saxild-Hansen and Per Chr. Hansen, May 4, 2010, DTU Compute.

% Reference: T. Elfving, T. Nikazad, and P. C. Hansen, Semi-convergence and
% relaxation parameters for a class of SIRT algorithms, Elec. Trans. on
% Numer. Anal., 37 (201), pp. 321-336.

if k(1) < 2
    error('The polynomial is not defined for k < 2')
end

% The starting guess for Newton's method.
z0 = 1;
% The number of Newton iterations.
kmax = 6;
z = zeros(length(k),1);

% Finds the root with Newton's method for all the values in k.
for i = 1:length(k)
    z(i) = myNewton(z0,kmax,k(i));
end

function [x,k] = myNewton(x0,kmax,grad)
% The Newton method specially defined for this problem.

% Initialize the parameters.
k = 0;
x = x0;
[f, df] = fung(x,grad);
h = f/df;

% Iterate using Newtoms method.
while  k < kmax
    k = k+1;
    x = x - h;
    [f, df] = fung(x,grad);
    h = f/df;
end

function [f, df] = fung(x,k)
% Uses Horner's algorithm to evaluate the ploynomial and its derivative.

C = [0 2*k-1 -ones(1,k-1)];
f = C(1);
df = 0;

for i = 2:length(C)
    df = df*x + f;
    f = f*x + C(i);
end
