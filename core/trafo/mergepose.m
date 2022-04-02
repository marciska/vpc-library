function g = mergepose(R,p) %#codegen
%MERGEPOSE returns homogeneuous representation from Rotation R and
%position p
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - R: Rotation matrix [size 3x3]
%   - p: Position vector [size 3]
%
% Outputs:
%   - g: Homogeneous representation / Pose matrix [size 4x4]
%
% Example commands:
%   g = mergepose(eye(3),[1 3 2]);
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2020
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% March 2021, Last modified: 2021-MAR-08
%
%------------- BEGIN CODE --------------

g = [R,       p(:); ...
     0, 0, 0,  1 ];

%------------- END CODE --------------
end