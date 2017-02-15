function [Mfun,sfun] = get_Mfun_sfun(sirt_method, A, m, n, M, w, s)

switch sirt_method
    
    case 'landweber'
        Mfun = @(XX) XX;
        sfun = @(XX) XX;
        
    case 'cimmino'
        sfun = @(XX) XX;
        % Define the M matrix.
        if isnan(M)
            % Calculate the norm of each row in A. This calculation can require a
            % lot of memory. The commented lines can be used instead; they are
            % slower, but use less memory!
            if ~isa(A,'function_handle')
                normAi = full(abs(sum(A.*A,2)));
                %       normAi = zeros(m,1);
                %       for i = 1:m
                %           ai = full(A(i,:));
                %           normAi(i) = norm(ai)^2;
                %       end
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
            %M = (1/m)*(1./normAi);
            
            I = (M == Inf);
            M(I) = 0;
        end
        Mfun = @(XX) M.*XX;
        
    case 'cav'
        sfun = @(XX) XX;
        % Define the M matrix.
        if isnan(M)
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
            
            % If the method is weigthed.
            if isnan(w)
                M = 1./normAs;
            else
                M = w./normAs;
            end
            I = (M == Inf);
            M(I) = 0;
        end
        Mfun = @(XX) M.*XX;
        
    case 'drop'
        % Define the M matrix.
        if isnan(M)
            
            % Calculate the norm of each row in A. This calculation can require
            % a lot of memory. The commented lines can be used instead; they
            % are slower, but use less memory!
            if ~isa(A,'function_handle')
                normAi = full(abs(sum(A.*A,2)));
                %           normAi = zeros(m,1);
                %           for i = 1:m
                %               ai = full(A(i,:));
                %               normAi(i) = norm(ai)^2;
                %           end
            else
                normAi = zeros(m,1);
                for i = 1:m
                    e = zeros(m,1);
                    e(i) = 1;
                    v = Afun(e,'transp');
                    normAi(i) = norm(v)^2;
                end
            end
            
            % If the method is weigthed.
            if isnan(w)
                M = 1./normAi;
            else
                M = w./normAi;
            end
            I = (M == Inf);
            M(I) = 0;
        end
        Mfun = @(XX) M.*XX;
        % Define the S matrix.
        if isnan(s)
            
            % Define the s vector and the M matrix.
            if ~isa(A,'function_handle')
                s = 1./sum(A~=0,1)';
            else
                s = zeros(n,1);
                for i=1:n
                    e = zeros(n,1); e(i) = 1;
                    Aj = Afun(e,'notransp');
                    s(i) = 1./sum(Aj~=0);
                end
            end
            I = (s == Inf);
            s(I) = 0;
        end
        sfun = @(XX) s.*XX;
        
    case 'sart'
        % Set diagonal of W = M if not given as input.
        if isnan(M)
            if ~isa(A,'function_handle')
                Aip = full(sum(abs(A),2));
            else
                Aip = abs(Afun(ones(n,1),'notransp'));
            end
            M = 1./Aip;
            I = (M == Inf);
            M(I) = 0;
        end
        Mfun = @(XX) M.*XX;
        % Set s-vector representing V = T if not given as input.
        if isnan(s)
            if ~isa(A,'function_handle')
                Apj = full(sum(abs(A),1))';
            else
                Apj = abs(Afun(ones(m,1),'transp'));
            end
            s = 1./Apj;
            I = (s == Inf);
            s(I) = 0;
        end
        sfun = @(XX) s.*XX;
        
        % Note on purpose weights not implemented for SART, as would modify
        % special property causing largest singular value to equal 1.
        
    otherwise
        error('SIRT method not defined.')
end