function dgwi = rbm(Vbwi,gwi) %#codegen
%RBM implements rigid body motion of a moving agent or target
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - Vbwi: Body velocity in world->object frame [size 6]
%   - gwi: state, homogeneous pose matrix g [size 4x4]
%
% Outputs:
%   - dgwi: derivative state, homogeneous pose matrix g [size 4x4]
%
% Example commands:
%   dgwi = rbm(ones(6,1),eye(4))
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

% Rigid Body Motion
dgwi = gwi * wedge(Vbwi);

%------------- END CODE --------------
end