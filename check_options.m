function options = check_options(options)

%% First check which fields are present. If not, set to default.

% if ~isfield(options,'ubound')
%     options.ubound = nan;
% end
% 
% if ~isfield(options,'nonneg')
%     options.nonneg = false;
% end

if ~isfield(options,'stoprule')
    options.stoprule.type = 'NO';
    options.stoprule.taudelta = nan;
else
    if ~isfield(options.stoprule,'type')
        options.stoprule.type = 'none';
    end
    if ~isfield(options.stoprule,'taudelta')
        options.stoprule.taudelta = nan;
    end
end

options.relaxpar
options.restart
options.restart.s1