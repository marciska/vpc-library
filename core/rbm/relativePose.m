function gBA = relativePose(gwA,gwB) %#codegen
%RELATIVEPOSE returns pose of agent A seen from agent B.
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - gwA: pose matrix of agent A seen from frame w [size 4x4]
%   - gwB: pose matrix of agent B seen from frame w [size 4x4]
%
% Outputs:
%   - gAB: pose matrix of agent B seen from frame B [size 4x4]
%
% Example commands:
%   gAB = relativePose(eye(4),eye(4))
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2021
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% March 2021, Last modified: 2021-MAR-09
%
%------------- BEGIN CODE --------------

% Relative pose
gBA = ginv(gwB)*gwA;

%-------------- END CODE ---------------
end