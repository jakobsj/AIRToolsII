function [k,K,rk,dk] = init_stoprules(stoprule,rxk,K,ncp_smooth)

% At this time, we know stoprule has been set. If no stoprule was given by
% user, now set to 'none'. If taudelta was missing, error has been thrown,
% so can now assume exists. We discard checks using isfield and ischar.
%
% Make sure to do check that field does exist previous to this, and that it
% is char.
%

% What to do about check that first iteration should be performed in DP?
% Not done in others, so different output arguments. Also call to return
% cannot terminate calling function, only this function. Probably move to
% calling function.
%                 % Check that the first iteration should be performed.
%                 rk = rxk;
%                 nrk = norm(rk);
%                 
%                 if nrk <= taudelta
%                     info = [2 k NaN];
%                     X = x0;
%                     return
%                 end % end the DP-rule.

% Make sure to check for taudelta
%                 if isfield(options.stoprule,'taudelta')
%                     taudelta = options.stoprule.taudelta;
%                 else
%                     error(['The factor taudelta must be specified when '...
%                         'using DP'])
%                 end

%                if isfield(options.stoprule,'taudelta')
%                     taudelta = options.stoprule.taudelta;
%                 else
%                     error(['The factor taudelta must be specified when '...
%                         'using ME'])
%                 end


% Defaults for k, dk and rk
k = 0;
rk = nan;
dk = nan;

switch upper(stoprule)
    case 'DP'
        % DP stopping rule: Nothing to do.
        
    case 'ME'
        % ME stopping rule.
        %rk = rxk;
        rk = nan(size(rxk));
        %K = K + 1;
        
    case 'NCP'
        % NCP stopping rule.
        dk = inf(ncp_smooth,1);
        K = [K max(K)+1];
        
    case 'NONE'
        % No stopping rule: Nothing to do.
        
    otherwise
        error('The chosen stopping rule is not valid.');
end