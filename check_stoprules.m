function [stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxpar, taudelta, k, kmax, rkm1, dk, res_dims)

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
                info.relaxpar = relaxpar;
            else
                info.stoprule = 0;
                info.finaliter = k;
                info.relaxpar = relaxpar;
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
            info.relaxpar = relaxpar;
        elseif k >= kmax
            stop = 1;
            info.stoprule = 0;
            info.finaliter = k;
            info.relaxpar = relaxpar;
        else
            rkm1 = rk;
        end % end the ME-rule.
        
    case 'NCP'
        % Normalized Cumulative Periodogram stopping rule.
        
        % Depending on residual being 1D or 2D, set up NCP
        switch length(res_dims)
            case 1
                % If 1D signal.
                m = length(rk);
                q = floor(m/2);
                rkh = fft(rk);
                pk = abs(rkh(1:q+1)).^2;
                c = cumsum(pk(2:end))/sum(pk(2:end));
                
            case 2
                % If 2D signal, reshape and do columnwise NCP -> mean.
                R = reshape(rk,res_dims);
                q = floor(res_dims(1)/2);
                RKH = fft(R);
                PK = abs(RKH(1:q+1,:)).^2;
                c = bsxfun(@rdivide,cumsum(PK(2:end,:)),sum(PK(2:end,:)));
                
            otherwise
                % Only implemented NCP for 1D and 2D signals.
                error(['NCP stopping rule only implemented for 1D and ',...
                    '2D residuals, the size of which should be given ',...
                    'using options.stoprule.res_dims.'])
        end
        
        % The cumulative white noise spectrum (spectrum is flat, hence
        % cumulative increasing linearly).
        c_white = (1:q)'./q;
        
        % Compute the norm of the difference between columns of c and
        % c_white. The below handles both 2D and 1D and is generalized from
        % the 1D statement   ncc = norm(c-c_white);
        ncc = mean(sqrt(sum(bsxfun(@minus,c,c_white).^2)));
        
        % Compare filtered new value with the previous filtered value.
        if  mean(dk) < mean([dk(2:end); ncc]) || k >= kmax
            stop = 1;
            if k ~= kmax
                info.stoprule = 1;
                info.finaliter = k;
                info.relaxpar = relaxpar;
            else
                info.stoprule = 0;
                info.finaliter = k;
                info.relaxpar = relaxpar;
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
            info.relaxpar = relaxpar;
        end % end stoprule type none
        
    otherwise
        error('Unknown stopping rule given.');
end % end stoprule type.