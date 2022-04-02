%% startup.m
% *Summary:* Loads necessary paths of the VPC library package
%
% -----------
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Supervisor:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2021
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% November 2021
%
% ------------- BEGIN CODE -------------

fprintf('Loading core VPC library...\n')

% add library root folder
libDir = fileparts(mfilename('fullpath'));
addpath(libDir)

% add paths to sub-dependencies
dirs = ["core", "templates"];
for d = dirs
    addpath(genpath(fullfile(libDir,d)))
end
clear d dirs libDir

% init gpml library
run('gpml/startup.m')

% clear workspace variables
clear d dirs libDir mydir

fprintf('[done]\n\n')

% -------------- END CODE --------------
