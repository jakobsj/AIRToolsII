function [stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, lambda, taudelta, k, kmax, rkm1, dk, res_dims)

% PCH: changed variable names for better consistency:
%      rk -> rkm1 (residual correcponding to iteration k-1)
%      rxk -> rk  (residual corresponding to iteration k)

% Only NCP uses dk: For all other stoprules, dk=nan is given as input and
% returned untouched as output.
% Only ME uses rkm1: For all other stoprules, rkm1=nan is given as input and
% returned untouched as output.

stop = 0;
info = nan;

switch upper(stoprule)
    case 'DP'
        % DP stopping rule.
        nrk = norm(rk);
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
        nrk = norm(rkm1);
        dME = rkm1'*(rkm1+rk)/2;  % PCH: still in accordance with redefined
                                  % ME stopping rule.
        
        %         if dME/nrk <= taudelta || k >= kmax
        %             stop = 1;
        %             if k ~= kmax
        %                 info = [3 k-1 lambda];
        %             else
        %                 info = [0 k-1 lambda];
        %             end
        % PCH replaced the above lines with the lines below.  NB: dette er
        %     en gammel kommentar fra en tidligere version.
        if dME/nrk <= taudelta
            stop = 1;
            info = [3 k lambda];  % PCH: changed "k-1" to "k" because I
                                  % redefined the ME stopping rule.
            % PCH: removed the warning below.
            % if k==1
            %     warning(['Initial vector satisfies stopping rule. '...
            %              'Returning initial vector in X.'])
            % end
        elseif k >= kmax
            stop = 1;
            info = [0 k lambda];
        else
            rkm1 = rk;
        end % end the ME-rule.
        
    case 'NCP'
        % NCP stopping rule.
        if length(res_dims) == 1
            m = length(rk);
            q = floor(m/2);
            c_white = (1:q)'./q;
            rkh = fft(rk);
            pk = abs(rkh(1:q+1)).^2;
%             c = zeros(q,1);
%             for index = 1:q
%                 c(index) = sum(pk(2:index+1))/sum(pk(2:end));
%             end
            % PCH: rewrote this using cumsum; perhaps move outside if-then-else.
            c = cumsum(pk(2:end))/sum(pk(2:end));
            
        else
            R = reshape(rk,res_dims);   % reshape(rxk,p,lt);
            RKH = fft2(R); % rkh = fft(r);
            q1 = floor(res_dims(1)/2); % PCH: a bit wasted effort to compute
            q2 = floor(res_dims(2)/2); % this stuff every time.
            PK = abs(RKH(1:q1+1,1:q2+1)).^2; % pk = abs(rkh(1:q+1)).^2; PCH: addded ".^2" !!!
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
%             c = zeros(q,1);
%             for index = 1:q
%                 c(index) = sum(pk(2:index+1))/sum(pk(2:end));
%             end
            % c_white = (1:q)'./q;  PCH: moved this line down.
            % PCH: rewrote this using cumsum; perhaps move outside if-then-else.
            c = cumsum(pk(2:end))/sum(pk(2:end));
        end
        
        c_white = (1:q)'./q;  % PCH: moved the above line to here.
                              % A bit wasted effort to compute it every time.
        
        ncc = norm(c-c_white); % PCH: avoid computing this norm twice.
        % if  mean(dk) < ncc || k >= kmax  % PCH: "norm(..." -> "ncc".
        % PCH: in the above implementation the unfiltered ncc is compared
        % with the previous filtered value; I think it is more correct to
        % compate a filtered new value with the previous filtered value:
        if  mean(dk) < mean([dk(2:end); ncc]) || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [1 k lambda];
      % PCH: changed "k-1" to "k" because I redefined the NCP stopping rule.
            else
                info = [0 k lambda];  % PCH: ditto.
            end
        else
            dk = [dk(2:end); ncc];  % PCH: "norm(..." -> "ncc".
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
