function [stop, info, rk, dk] = check_stoprules(...
    stoprule, rxk, lambda, taudelta, k, kmax, rk, dk)

% Only NCP uses dk: For all other stoprules, dk=nan is given as input and
% returned untouched as output.
% Only ME uses rk: For all other stoprules, rk=nan is given as input and
% returned untouched as output.

switch upper(stoprule)
    case 'DP'
        % DP stopping rule.
        nrk = norm(rxk);
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [2 k lambda];
            else
                info = [0 k lambda];
            end
        end % end the DP-rule.
    
    case 'ME'
        % ME stopping rule.
        nrk = norm(rk);
        dME = rk'*(rk+rxk)/2;
        
        if dME/nrk <= taudelta || k >= kmax
            stop = 1;
            
            if k ~= kmax
                info = [3 k-1 lambda];
            else
                info = [0 k-1 lambda];
            end
        else
            rk = rxk;
        end % end the ME-rule.
        
    case 'NCP'
        % NCP stopping rule.
        m = length(rxk);
        q = floor(m/2);
        c_white = (1:q)'./q;
        rkh = fft(rxk);
        pk = abs(rkh(1:q+1)).^2;
        c = zeros(q,1);
        for index = 1:q
            c(index) = sum(pk(2:index+1))/sum(pk(2:end));
        end
        
        if dk < norm(c-c_white) || k >= kmax
            stop = 1;            
            if k ~= kmax
                info = [1 k-1 lambda];
            else
                info = [0 k-1 lambda];
            end            
        else
            dk = norm(c-c_white);
        end % end NCP-rule.
        
    case 'NONE'
        % No stopping rule.
        if k >= kmax 
            stop = 1;
            info = [0 k lambda];
        end % end stoprule type none
        
    otherwise
        error('Unknown stopping rule given.');
end % end stoprule type.
        
        
        
% Stopping rules:
%     if strncmpi(stoprule,'DP',2)
%         % DP stopping rule.
%         nrk = norm(rxk);
%         
%         if nrk <= taudelta || k >= kmax
%             stop = 1;
%             if k ~= kmax
%                 info = [2 k lambda];
%             else
%                 info = [0 k lambda];
%             end
%         end % end the DP-rule.
%         
%     elseif strncmpi(stoprule,'ME',2)
%         % ME stopping rule.
%         nrk = norm(rk);
%         dME = rk'*(rk+rxk)/2;
%         
%         if dME/nrk <= taudelta || k >= kmax
%             stop = 1;
%             
%             if k ~= kmax
%                 info = [3 k-1 lambda];
%             else
%                 info = [0 k-1 lambda];
%             end
%         else
%             rk = rxk;
%         end % end the ME-rule.
%         
%     elseif strncmpi(stoprule,'NC',2)
%         % NCP stopping rule.
%         m = length(rxk);
%         q = floor(m/2);
%         c_white = (1:q)'./q;
%         rkh = fft(rxk);
%         pk = abs(rkh(1:q+1)).^2;
%         c = zeros(q,1);
%         for index = 1:q
%             c(index) = sum(pk(2:index+1))/sum(pk(2:end));
%         end
%         
%         if dk < norm(c-c_white) || k >= kmax
%             stop = 1;            
%             if k ~= kmax
%                 info = [1 k-1 lambda];
%             else
%                 info = [0 k-1 lambda];
%             end            
%         else
%             dk = norm(c-c_white);
%         end % end NCP-rule.
%         
%     elseif strncmpi(stoprule,'NO',2)
%         % No stopping rule.
%         if k >= kmax 
%             stop = 1;
%             info = [0 k lambda];
%         end
%     end % end stoprule type.
