function logRot = logR(Rot) %#codegen
%LOGR is the inverse function to exponential rotation matrix
%exp{xi^{wedge}*theta}
%
% Detailed Explanation:
%   Returns the argument of exp{xi^{wedge}*theta}.
%   In other words, the matrix xi^{wedge}*theta is returned.
%
% -----------
%
% Inputs:
%   - Rot: rotation matrix [size 3x3]
%
% Outputs:
%   - logRot: matrix [size 3x3]
%
% Example commands:
%   logRot = logR(eye(3));
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita Lab, University of Tokyo, 2021
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% March 2021, Last modified: 2021-MAR-09
%
%------------- BEGIN CODE --------------

if all(Rot == eye(3), 'all') % Special case: Rot is identity matrix
    logRot = zeros(3); % No rotation
else % Normal case
    cosR_ = .5*(trace(Rot)-1);
    cosR = max(-1, min(cosR_, 1)); % -1 ≤ cosR ≤ 1
    phiRot = acos(cosR); % angle of rotation
    logRot = phiRot/(2*sin(phiRot))*(Rot-Rot');
end

%-------------- END CODE ---------------
end