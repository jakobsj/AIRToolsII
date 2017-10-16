function [stop, info, rkm1, dk] = check_stoprules(...
    stoprule, rk, relaxpar, taudelta, k, kmax, rkm1, dk, res_dims)
%CHECK_STOPRULES  Check if stopping rule criteria are met
%
%   [stop, info, rkm1, dk] = check_stoprules(...
%    stoprule, rk, relaxpar, taudelta, k, kmax, rkm1, dk, res_dims)
%
% From information available at current iteration of an ART, CART or SIRT
% method (and for some stopping rules previous iterates) determine if the
% chosen stopping rule is satisfied in order to abort iteration.
% 
% Input:
%    stoprule     Name of stopping rule: 'none', 'DP', 'ME', or 'NCP'.
%    rk           Current residual vector.
%    relaxpar     Current relaxation parameter.
%    taudelta     Parameter for use in DP and ME stopping rules.
%    k            The current iteration number.
%    kmax         The maximal number of iterations to run.
%    rkm1         Residual vector from previous iteration. Used only by ME.
%    dk           Vector with values to filter/average over from previous 
%                 iterations. Used only by NCP.
%    res_dims     Dimensions of residual. Used only by NCP.
%    
% Output:
%    stop         1 if stopping rule is met, otherwise 0.
%    info         Struct with fields 
%       stoprule = 0 : stopped by maximum number of iterations
%                  1 : stopped by NCP-rule
%                  2 : stopped by DP-rule
%                  3 : stopped by ME-rule.
%       finaliter    : no. of iterations in total.
%       relaxpar     : the chosen relaxation parameter.
%    rkm1         Current residual saved for use as previous in next iter.
%    dk           Filter vector updated by replacing oldest element by
%                 the one from current iteration.
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

% Defaults.
stop = 0;
info = struct;

% Split by stopping rule chosen.
switch upper(stoprule)
    
    case 'DP'
        % Discrepancy Principle stopping rule.
        nrk = norm(rk);
        if nrk <= taudelta || k >= kmax
            stop = 1;
            if k ~= kmax
                % Discrepancy principle stopping rule satisfied.
                info.stoprule = 2;
                info.finaliter = k;
                info.relaxpar = relaxpar;
            else
                % Maximum number of iterations reached.
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
            % Monotone Error stopping rule satisfied.
            stop = 1;
            info.stoprule = 3;
            info.finaliter = k;
            info.relaxpar = relaxpar;
        elseif k >= kmax
            % Maximum number of iterations reached.
            stop = 1;
            info.stoprule = 0;
            info.finaliter = k;
            info.relaxpar = relaxpar;
        else
            % No stopping rule satisfied, prepare for next iteration.
            rkm1 = rk;
        end % end the ME-rule.
        
    case 'NCP'
        % Normalized Cumulative Periodogram stopping rule.
        
        % Depending on residual being 1D or 2D, set up NCP.
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
        % c_white. The below handles both 2D and 1D and is generalized
        % from the 1D statement:  ncc = norm(c-c_white).
        ncc = mean(sqrt(sum(bsxfun(@minus,c,c_white).^2)));
        
        % Compare filtered new value with the previous filtered value.
        if  max(dk) < max([dk(2:end); ncc]) || k >= kmax
            stop = 1;
            if k ~= kmax
                % NCP stopping rule satisfied.
                info.stoprule = 1;
                info.finaliter = k;
                info.relaxpar = relaxpar;
            else
                % Maximum number of iterations reached.
                info.stoprule = 0;
                info.finaliter = k;
                info.relaxpar = relaxpar;
            end
        else
            % No stopping rule satisfied, prepare for next iteration:
            % Replace oldest value by shifting all and storing newest.
            dk = [dk(2:end); ncc];
        end % end NCP-rule.
        
    case 'NONE'
        % No stopping rule.
        if k >= kmax
            % Maximum number of iterations reached.
            stop = 1;
            info.stoprule = 0;
            info.finaliter = k;
            info.relaxpar = relaxpar;
        end % end stoprule type none.
        
    otherwise
        error('Unknown stopping rule given.');
end % end stoprule type.