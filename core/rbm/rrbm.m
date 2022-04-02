function dgco = rrbm(Vbwc,Vbwo,gco) %#codegen
%RRBM returns the relative rigid body motion of a target seen from a Camera
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - Vbwo: Body velocity in world->Target frame [size 6]
%   - Vbwc: Body velocity in world->Camera frame [size 6]
%   - gco: state, homogeneous pose matrix g [size 4x4]
%
% Outputs:
%   - dgco: derivative state, homogeneous pose matrix g [size 4x4]
%
% Example commands:
%   dgco = rrbm(ones(6,1),ones(6,1),eye(4))
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

% Relative Rigid Body Motion
dgco = -wedge(Vbwc)*gco + gco*wedge(Vbwo);

%------------- END CODE --------------
end