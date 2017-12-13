function AIRToolsII_setup(choice)
%AIRToolsII_setup  Function that sets up search paths to AIR Tools II
%
% AIRToolsII_setup(choice)
%
% Run this function from the top directory of AIR Tools II, to include the
% package in MATLAB's path.  Use the input choice in this way:
%   choice = 'temporary': the path is temporarily updated and used for the
%            current MATLAB session only - this is the default choice,
%   choice = 'permanent': the path is updated and saved for both the
%            current and all future MATLAB sessions.
% In the latter case, if you want to remove the paths associated with
% AIR Tools II, use rmpath or use the "Set Path" tool in the "Home" tab
% of the main MATLAB window.

% Code written by: Per Christian Hansen, Jakob Sauer Jorgensen, and 
% Maria Saxild-Hansen, 2010-2017.

% This file is part of the AIR Tools II package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2017 Per Christian Hansen, Technical University of Denmark and
% Jakob Sauer Jorgensen, University of Manchester.

if nargin==0, choice = 'temporary'; end

addpath(genpath(fileparts(mfilename('fullpath'))));

if strcmp(choice,'permanent')
    status = savepath;
    if status==1
        warning('AIR Tools II was added to the MATLAB search path for the current session only. Adding it permanently failed, probably due to a write permission issue. It is possible to manually add AIR Tools II permanently to your search path, but may require consulting your system administrator. Alternatively, you can re-run AIRToolsII_setup without an input parameter in each new MATLAB session where you want to use AIR Tools II.')
    end
end