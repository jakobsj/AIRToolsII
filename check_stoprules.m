function [stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, lambda, taudelta, k, kmax, rkm1, dk, res_dims)

% Only NCP uses dk: For all other stoprules, dk=nan is given as input 
% and % returned untouched as output.
% Only ME uses rkm1: For all other stoprules, rkm1=nan is given as input 
% and returned untouched as output.

stop = 0;

info = struct;

switch upper(stoprule)
    
    case 'DP'
        % Discrepancy Principle stopping rule.
        nrk = norm(rk);
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k ~= kmax
                info.stoprule = 2;
                info.finaliter = k;
                info.lambda = lambda;
            else
                info.stoprule = 0;
                info.finaliter = k;
                info.lambda = lambda;
            end
        end % end the DP-rule.
    
    case 'ME'
        % Monotone Error stopping rule.
        nrk = norm(rkm1);
        dME = rkm1'*(rkm1+rk)/2;

        if dME/nrk <= taudelta
            stop = 1;
            info.stoprule = 3;
            info.finaliter = k;
            info.lambda = lambda;
        elseif k >= kmax
            stop = 1;
            info.stoprule = 0;
            info.finaliter = k;
            info.lambda = lambda;
        else
            rkm1 = rk;
        end % end the ME-rule.
        
    case 'NCP'
        % Normalized Cumulative Periodogram stopping rule.
        if length(res_dims) == 1
            m = length(rk);
            q = floor(m/2);
            rkh = fft(rk);
            pk = abs(rkh(1:q+1)).^2;
            % PCH: rewrote this using cumsum; 
            % perhaps move outside if-then-else.
            c = cumsum(pk(2:end))/sum(pk(2:end));
        else
            R = reshape(rk,res_dims);   % reshape(rxk,p,lt);
            RKH = fft2(R); % rkh = fft(r);
            q1 = floor(res_dims(1)/2); % PCH: a bit wasted effort to 
            q2 = floor(res_dims(2)/2); % compute this stuff every time.
            PK = abs(RKH(1:q1+1,1:q2+1)).^2; % pk = abs(rkh(1:q+1)).^2; 
            P = zeros(size(PK));  % PCH: note to self, this is faster than
                                  % [I,J] = meshgrid(0:n-1); P = I.^2+J.^2;
                                  % and uses less memory.
            for i=1:size(P,1)
                for j=1:size(P,2)
                    P(i,j) = (i-1)^2+(j-1)^2;
                end
            end
            [~,I] = sort(P(:));
            pk = PK(I);
            q = length(pk)-1;
            % PCH: rewrote this using cumsum; 
            % perhaps move outside if-then-else.
            c = cumsum(pk(2:end))/sum(pk(2:end));
        end
        
        c_white = (1:q)'./q;  % PCH: moved the above line to here.
                              % A bit wasted effort to compute every time.
        
        ncc = norm(c-c_white); % PCH: avoid computing this norm twice.
        
        % Compare filtered new value with the previous filtered value.
        if  mean(dk) < mean([dk(2:end); ncc]) || k >= kmax
            stop = 1;
            if k ~= kmax
                info.stoprule = 1;
                info.finaliter = k;
                info.lambda = lambda;
            else
                info.stoprule = 0;
                info.finaliter = k;
                info.lambda = lambda;
            end
        else
            % Replace oldest value by shifting all and storing newest.
            dk = [dk(2:end); ncc];
        end % end NCP-rule.
        
    case 'NONE'
        % No stopping rule.
        if k >= kmax 
            stop = 1;
            info.stoprule = 0;
            info.finaliter = k;
            info.lambda = lambda;
        end % end stoprule type none
        
    otherwise
        error('Unknown stopping rule given.');
end % end stoprule type.