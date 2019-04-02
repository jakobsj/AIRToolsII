function [Mfun,Dfun,Mflag,Dflag] = get_mfun_dfun(sirt_method, A, m, n)
%GET_MFUN_DFUN  Computes matrices M and D for SIRT methods
%
%   [Mfun,Dfun,Mflag,Dflag] = get_mfun_dfun(sirt_method, A, m, n)
%
% Set up the matrices M and D characterizing a SIRT method.
% 
% Input:
%    sirt_method  Name of SIRT method: 'landweber', 'cimmino', 'cav', 
%                 'drop', or 'sart'.
%    A            The forward operator, either as matrix or function handle.
%    m            Number of rows in A.
%    n            Number of columns in A.
%    
% Output:
%    Mfun         Function handle which applies the matrix-vector
%                 multiplication of M on input, i.e., Mfun(v) is M*v.
%    Dfun         Function handle which applies the matrix-vector
%                 multiplication of D on input, i.e., Dfun(v) is D*v.
%    Mflag        Integer stating how the matrix operation within the Mfun
%                 function handle is represented. 0: the matrix M is the
%                 identity represented as simply returning the input. 1:
%                 a vector of weights corresponding to the diagonal of M.
%                 2: a matrix.
%    Dflag        Integer stating how the matrix operation within the Dfun
%                 function handle is represented. 0: the matrix D is the
%                 identity represented as simply returning the input. 1:
%                 a vector of weights corresponding to the diagonal of D.
%                 2: a matrix.
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

% Set block size for row norm computations; adjust if necessary.
B = 200;

% Split between built-in and custom methods depending on type of sirt_method.
if ischar(sirt_method)
    switch sirt_method
        
        case 'landweber'
            
            % Both M and D are identity, so do nothing.
            Mfun = @(XX) XX;
            Dfun = @(XX) XX;
            
        case 'cimmino'
            
            % D is identity, so do nothing.
            Dfun = @(XX) XX;
            
            % Define the M matrix. Compute the norm of each row in A.
            normAi = zeros(m,1);
            if isnumeric(A)
                for J = 1:ceil(m/B)
                    I = 1+(J-1)*B : min(m,J*B);
                    normAi(I) = sum(A(I,:).*A(I,:),2);
                end
            elseif isa(A,'function_handle')
                for i = 1:m
                    e = zeros(m,1); e(i) = 1;
                    v = A(e,'transp');
                    normAi(i) = norm(v)^2;
                end
            else
                for i = 1:m
                    e = zeros(m,1); e(i) = 1;
                    v = A'*e;
                    normAi(i) = norm(v)^2;
                end
            end
            
            % Generate diagonal entries and fix division by zero.
            M = 1/m*(1./normAi);
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
        case 'cav'
            
            % D is identity, so do nothing.
            Dfun = @(XX) XX;
            
            % Define the M matrix.
            s = zeros(n,1);
            normAs = zeros(m,1);
            if isnumeric(A)
                for I = 1:ceil(n/B)
                    J = 1+(I-1)*B : min(n,I*B);
                    s(J) = sum(A(:,J)~=0,1)';
                end
                s = spdiags(s,0,n,n);
                for J = 1:ceil(m/B)
                    I = 1+(J-1)*B : min(m,J*B);
                    normAs(I) = sum((A(I,:).*A(I,:))*s,2);
                end
            elseif isa(A,'function_handle')
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = A(e,'notransp');
                    s(i) = 1/sum(Aj~=0);
                end
                for i=1:m
                    e = zeros(m,1); e(i) = 1;
                    Ai = A(e,'transp');
                    normAs(i) = (s'*Ai.^2);
                end
            else
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = A*e;
                    s(i) = 1/sum(Aj~=0);
                end
                for i=1:m
                    e = zeros(m,1); e(i) = 1;
                    Ai = A'*e;
                    normAs(i) = (s'*Ai.^2);
                end
            end
            
            % Generate diagonal entries and fix division by zero.
            M = 1./normAs;
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
        case 'drop'
            
            % Define the M matrix; same as in Cimmino.
            % Compute the norm of each row in A.
            normAi = zeros(m,1);
            if isnumeric(A)
                for J = 1:ceil(m/B)
                    I = 1+(J-1)*B : min(m,J*B);
                    normAi(I) = sum(A(I,:).*A(I,:),2);
                end
            elseif isa(A,'function_handle')
                for i = 1:m
                    e = zeros(m,1); e(i) = 1;
                    v = A(e,'transp');
                    normAi(i) = norm(v)^2;
                end
            else
                for i = 1:m
                    e = zeros(m,1); e(i) = 1;
                    v = A'*e;
                    normAi(i) = norm(v)^2;
                end
            end
            
            % Generate diagonal entries and fix division by zero.
            M = 1./normAi;
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
            % Define diagonal s of the D matrix.
            s = zeros(n,1);
            if isnumeric(A)
                for I = 1:ceil(n/B)
                    J = 1+(I-1)*B : min(n,I*B);
                    s(J) = 1./sum(A(:,J)~=0,1)';
                end
            elseif isa(A,'function_handle')
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = A(e,'notransp');
                    s(i) = 1./sum(Aj~=0);
                end
            else
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = A*e;
                    s(i) = 1./sum(Aj~=0);
                end
            end
            
            % Fix divisions by zero.
            I = (s == Inf);
            s(I) = 0;
            
            % Store D in function handle.
            Dfun = @(XX) s.*XX;
            
        case 'sart'
            % Compute diagonal of M.
            if isnumeric(A)
                Aip = zeros(m,1);
                for J = 1:ceil(m/B)
                    I = 1+(J-1)*B : min(m,J*B);
                    Aip(I) = sum(abs(A(I,:)),2);
                end
            elseif isa(A,'function_handle')
                Aip = abs(A(ones(n,1),'notransp'));
            else
                Aip = abs(A*ones(n,1));
            end
            M = 1./Aip;
            
            % Fix divisions by zero.
            I = (M == Inf);
            M(I) = 0;
            
            % Store M in function handle.
            Mfun = @(XX) M.*XX;
            
            %  Compute diagonal of D.
            if isnumeric(A)
                Apj = zeros(n,1);
                for I = 1:ceil(n/B)
                    J = 1+(I-1)*B : min(n,I*B);
                    Apj(J) = sum(abs(A(:,J)),1);
                end
            elseif isa(A,'function_handle')
                Apj = abs(A(ones(m,1),'transp'));
            else
                Apj = abs(A'*ones(m,1));
            end
            D = 1./Apj;
            
            % Fix divisions by zero.
            I = (D == Inf);
            D(I) = 0;
            
            % Store D as function handle.
            Dfun = @(XX) D.*XX;
            
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
