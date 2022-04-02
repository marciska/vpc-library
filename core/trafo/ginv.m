function g_inv = ginv(g) %#codegen
%GINV returns inverse of pose matrix g
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - g: pose matrix [size 4x4]
%
% Outputs:
%   - g_inv: inverse of pose matrix [size 4x4]
%
% Example commands:
%   g_inv = ginv(mergepose(eye(3),[1 3 2]));
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

% Extract position & orientation from g
[R,p] = splitpose(g);

% Invert pose
g_inv = [    R',     -R'*p(:); ...
         zeros(1,3),     1   ];

%-------------- END CODE ---------------
end