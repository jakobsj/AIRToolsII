function [Mfun,Tfun] = get_Mfun_Tfun(sirt_method, A, m, n, w, savememory)

% New option "savememory" to allow user to choose whether to use the
% standard fast-but-memory-demanding computation of M and T, or the
% slower-but-memory-efficient approach, more similar to the matrix-free.

switch sirt_method
    
    case 'landweber'
        
        % Both M and T are identity, so do nothing.
        Mfun = @(XX) XX;
        Tfun = @(XX) XX;
        
    case 'cimmino'
        
        % T is identity, so do nothing.
        Tfun = @(XX) XX;
        
        % Define the M matrix.
        % Calculate the norm of each row in A. 
        % Default (savememory=false) works on whole matrix at once, which
        % is fast but memory-demanding.
        % Alternative is savememory=true which works row-by-row, same as
        % in the function_handle/matrix-free approach.
        if ~isa(A,'function_handle')
            if savememory
                normAi = zeros(m,1);
                for i = 1:m
                    ai = full(A(i,:));
                    normAi(i) = norm(ai)^2;
                end
            else
                normAi = full(abs(sum(A.*A,2)));
            end
        else
            normAi = zeros(m,1);
            for i = 1:m
                e = zeros(m,1);
                e(i) = 1;
                v = A(e,'transp');
                normAi(i) = norm(v)^2;
            end
        end
        
        % If the method is weigthed.
        if isnan(w)
            M = 1/m*(1./normAi);
        else
            M = 1/m*(w./normAi);
        end

        % Fix any divisions be zero.
        I = (M == Inf);
        M(I) = 0;
        
        % Store M in function handle.
        Mfun = @(XX) M.*XX;
        
    case 'cav'
        
        % T is identity, so do nothing.
        Tfun = @(XX) XX;
        
        % Define the M matrix.
        if ~isa(A,'function_handle')
            if savememory
                error('not implemented')
            else
                s = sum(A~=0,1)';
                s = spdiags(s,0,n,n);
                normAs = full(sum(A.^2*s,2));
            end
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
        
        % If the method is weighted.
        if isnan(w)
            M = 1./normAs;
        else
            M = w./normAs;
        end
        
        % Fix divisions by zero.
        I = (M == Inf);
        M(I) = 0;
        
        % Store M in function handle.
        Mfun = @(XX) M.*XX;
        
    case 'drop'
        
        % Define the M matrix. Same as in Cimmino.
        % Calculate the norm of each row in A. 
        % Default (savememory=false) works on whole matrix at once, which
        % is fast but memory-demanding.
        % Alternative is savememory=true which works row-by-row, same as
        % in the function_handle/matrix-free approach.
        if ~isa(A,'function_handle')
            if savememory
                normAi = zeros(m,1);
                for i = 1:m
                    ai = full(A(i,:));
                    normAi(i) = norm(ai)^2;
                end
            else
                normAi = full(abs(sum(A.*A,2)));
            end
        else
            normAi = zeros(m,1);
            for i = 1:m
                e = zeros(m,1);
                e(i) = 1;
                v = Afun(e,'transp');
                normAi(i) = norm(v)^2;
            end
        end
        
        % If the method is weighted.
        if isnan(w)
            M = 1./normAi;
        else
            M = w./normAi;
        end
        
        % Fix divisions by zero.
        I = (M == Inf);
        M(I) = 0;
        
        % Store M in function handle.
        Mfun = @(XX) M.*XX;
        
        % Define the T matrix.
        % Define the s vector and the M matrix.
        if ~isa(A,'function_handle')
            if savememory
                error('not implemented')
            else
                s = 1./sum(A~=0,1)';
            end
        else
            s = zeros(n,1);
            for i=1:n
                e = zeros(n,1); e(i) = 1;
                Aj = Afun(e,'notransp');
                s(i) = 1./sum(Aj~=0);
            end
        end
        
        % Fix divisions by zero.
        I = (s == Inf);
        s(I) = 0;
        
        % Store T in function handle.
        Tfun = @(XX) s.*XX;
        
    case 'sart'
        % Set diagonal of W = M if not given as input. W is the notation of
        % the original SART article.
        if ~isa(A,'function_handle')
            if savememory
                error('not implemented')
            else
                Aip = full(sum(abs(A),2));
            end
        else
            Aip = abs(Afun(ones(n,1),'notransp'));
        end
        M = 1./Aip;
        
        % Fix divisions by zero.
        I = (M == Inf);
        M(I) = 0;
        
        % Store M in function handle.
        Mfun = @(XX) M.*XX;
        
        % Set s-vector representing V = T if not given as input. V is the
        % notation from the original SART article.
        if ~isa(A,'function_handle')
            if savememory
                error('not implemented')
            else
                Apj = full(sum(abs(A),1))';
            end
        else
            Apj = abs(Afun(ones(m,1),'transp'));
        end
        s = 1./Apj;
        
        % Fix divisions by zero.
        I = (s == Inf);
        s(I) = 0;
        
        % Store T as function handle.
        Tfun = @(XX) s.*XX;
        
        % Note on purpose weights not implemented for SART, as would modify
        % special property causing largest singular value to equal 1.
        
    otherwise
        error('SIRT method not defined.')
end