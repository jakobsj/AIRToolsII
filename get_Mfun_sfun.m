function [Mfun,sfun] = get_Mfun_sfun(sirt_method, A, m, M, w)

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
                
        
    otherwise
        error('SIRT method not defined.')
end