function [Mfun,Dfun,Mflag,Dflag] = get_mfun_dfun(sirt_method, A, m, n)
%GET_MFUN_DFUN Aux. function to set up M and D matrices for SIRT methods
%
%   [Mfun,Dfun,Mflag,Dflag] = get_mfun_dfun(sirt_method, A, m, n)
%
% Set up the matrices M and D characterizing a SIRT method.
% 
% Input:
%    sirt_method  Name of SIRT method: 'landweber', 'cimmino', 'cav', 
%                 'drop', or 'sart'.
%    A            The forward operator A, either as matrix or function.
%    m            Number of rows in A.
%    n            Number of columns in A.
%    
% Output:
%    Mfun         Function handle which applies the matrix-vector
%                 multiplication of M on input Mfun(v) is M*v.
%    Dfun         Function handle which applies the matrix-vector
%                 multiplication of D on input Dfun(v) is D*v.
%    Mflag        Integer stating how the matrix operation within the Mfun
%                 function handle is represented. 0: The matrix M is the
%                 identity represented as simply returning the input. 1:
%                 A vector of weights corresponding to the diagonal of M.
%                 2: A matrix.
%    Dflag        Integer stating how the matrix operation within the Dfun
%                 function handle is represented. 0: The matrix D is the
%                 identity represented as simply returning the input. 1:
%                 A vector of weights corresponding to the diagonal of D.
%                 2: A matrix.
%
% See also: art.m, cart.m, sirt.m

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, DTU Compute, 2010-2017.

% This file is part of the AIR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package. 
% 
% Copyright 2017 Per Christian Hansen & Jakob Sauer Jorgensen, DTU Compute


% Split between builtin and custom methods depending on type of sirt_method
if ischar(sirt_method)
    switch sirt_method
        
        case 'landweber'
            
            % Both M and D are identity, so do nothing.
            Mfun = @(XX) XX;
            Dfun = @(XX) XX;
            
        case 'cimmino'
            
            % D is identity, so do nothing.
            Dfun = @(XX) XX;
            
            % Define the M matrix.
            % Calculate the norm of each row in A.
            if ~isa(A,'function_handle')
                normAi = full(abs(sum(A.*A,2)));
            else
                normAi = zeros(m,1);
                for i = 1:m
                    e = zeros(m,1);
                    e(i) = 1;
                    v = A(e,'transp');
                    normAi(i) = norm(v)^2;
                end
            end
            
            % Generate diagonal entries.
            M = 1/m*(1./normAi);
            
            % Fix any divisions be zero.
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
        case 'cav'
            
            % D is identity, so do nothing.
            Dfun = @(XX) XX;
            
            % Define the M matrix.
            if ~isa(A,'function_handle')
                s = sum(A~=0,1)';
                s = spdiags(s,0,n,n);
                normAs = full(sum(A.^2*s,2));
            else
                s = zeros(n,1);
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = A(e,'notransp');
                    s(i) = 1/sum(Aj~=0);
                end
                normAs = zeros(m,1);
                for i=1:m
                    e = zeros(m,1); e(i) = 1;
                    Ai = A(e,'transp');
                    normAs(i) = (s'*Ai.^2);
                end
            end
            
            % Generate diagonal entries.
            M = 1./normAs;
            
            % Fix divisions by zero.
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
        case 'drop'
            
            % Define the M matrix. Same as in Cimmino.
            % Calculate the norm of each row in A.
            if ~isa(A,'function_handle')
                normAi = full(abs(sum(A.*A,2)));
            else
                normAi = zeros(m,1);
                for i = 1:m
                    e = zeros(m,1);
                    e(i) = 1;
                    v = A(e,'transp');
                    normAi(i) = norm(v)^2;
                end
            end
            
            % Generate diagonal entries.
            M = 1./normAi;
            
            % Fix divisions by zero.
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
            % Define the D matrix.
            % Define the s vector and the M matrix.
            if ~isa(A,'function_handle')
                s = 1./sum(A~=0,1)';
            else
                s = zeros(n,1);
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = A(e,'notransp');
                    s(i) = 1./sum(Aj~=0);
                end
            end
            
            % Fix divisions by zero.
            I = (s == Inf);
            s(I) = 0;
            
            % Store D in function handle.
            Dfun = @(XX) s.*XX;
            
        case 'sart'
            % Set diagonal of W = M if not given as input. W is the notation of
            % the original SART article.
            if ~isa(A,'function_handle')
                Aip = full(sum(abs(A),2));
            else
                Aip = abs(A(ones(n,1),'notransp'));
            end
            
            % Generate diagonal entries.
            M = 1./Aip;
            
            % Fix divisions by zero.
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
            % Set s-vector representing V = D if not given as input. V is the
            % notation from the original SART article.
            if ~isa(A,'function_handle')
                Apj = full(sum(abs(A),1))';
            else
                Apj = abs(A(ones(m,1),'transp'));
            end
            s = 1./Apj;
            
            % Fix divisions by zero.
            I = (s == Inf);
            s(I) = 0;
            
            % Store D as function handle.
            Dfun = @(XX) s.*XX;
            
        otherwise
            error('SIRT method not defined.')
    end
    
    % Set Mflag and Dflag to 1 to show that for all built-in SIRT methods,
    % M and D are represented as a vector.
    Mflag = 1;
    Dflag = 1;
    
elseif isstruct(sirt_method) % Then it is a custom method.
    % Possible to pass in custom SIRT method given by struct with fields M
    % and D holding matrices or vectors (of diagonal elements instead of
    % string input.
    if isfield(sirt_method,'M') && ~isempty(sirt_method.M)
        if isnumeric(sirt_method.M)
            % Now we know that field M exists, is numeric and not empty.
            % Decide if vector.
            if isvector(sirt_method.M)
                Mfun = @(in) sirt_method.M.*in;
                Mflag = 1;
            else
                % If the matrix is diagonal, extract diagonal for
                % efficiency, otherwise keep matrix.
                if isdiag(sirt_method.M)
                    Mfun = @(in) diag(sirt_method.M).*in;
                    Mflag = 1;
                else
                    Mfun = @(in) sirt_method.M*in;
                    Mflag = 2;
                end
            end
        else
            error('M must be given as an array, either of dimension 0, 1 or 2.')
        end
    else % If no M, use identity.
        Mfun = @(in) in;
        Mflag = 0;
    end
    
    % Same for D.
    if isfield(sirt_method,'D') && ~isempty(sirt_method.D)
        if isnumeric(sirt_method.D)
            % Now we know that field D exists, is numeric and not empty.
            % Decide if vector.
            if isvector(sirt_method.D)
                Dfun = @(in) sirt_method.D.*in;
                Dflag = 1;
            else
                % If the matrix is diagonal, extract diagonal for
                % efficiency, otherwise keep matrix.
                if isdiag(sirt_method.D)
                    Dfun = @(in) diag(sirt_method.D).*in;
                    Dflag = 1;
                else
                    Dfun = @(in) sirt_method.D*in;
                    Dflag = 2;
                end
            end
        else
            error('D must be given as an array, either of dimension 0, 1 or 2.')
        end
    else % If no D, use identity.
        Dfun = @(in) in;
        Dflag = 0;
    end
else
    error(['First input must be a string with the SIRT method name or ',...
        'a struct with fields M and D characterizing a custom method.']);
end
