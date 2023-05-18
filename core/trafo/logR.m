function xiWtheta = logR(R) %#codegen
%%LOGR returns argument of exponential rotation matrix exp{xi^{∧}*θ}.
%%In other words, the matrix xi^{∧}*θ is returned.
%
% Syntax:   XIWTHETA = LOGR(R)
%
% Inputs:
%
%   - R: rotation matrix [size 3x3]
%
% Outputs:
%
%   - XIWTHETA: matrix [size 3x3]
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, March 2021
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

%% new implementation
% BUG: the section below has an issue when rot == pi. In practice, this
% will rarily create an issue because the numerics is never exactly pi in
% our simulation. Nevertheless, we should switch to 'rotm2axang' since it
% considers this corner case. Problem is, that rotm2axang creates some
% simulink errors about variable signals... we'll solve this at a later
% version of the VPC library. Thus, we'll comment it out for now here.

% axang = rotm2axang(R); % axang = [xi, theta]
% xi = axang(1:3); theta = axang(4);
% xiWtheta = wedge(xi*theta);


%% old, deprecated implementation
% BUG: Check note above. Below code has an error when rot == pi

if all(R == eye(3), 'all') % Special case: Rot is identity matrix
    xiWtheta = zeros(3); % No rotation
else % Normal case
    cosR_ = .5*(trace(R)-1);
    cosR = max(-1, min(cosR_, 1)); % -1 ≤ cosR ≤ 1
    phiR = acos(cosR); % angle of rotation
    % TODO: There is a problem here. If the rotation is pi, then R-R'
    % is zero, and thus logRot is completely zeros. This falsely shows
    % indicates that the rotation angle is zero, since the axis of rotation
    % cannot be determined. Is there a way to determine the rotation axis?
    dR = R-R';
%     if all(dR == 0, 'all')
%         % handle special case here since rotation axis cannot be determined
%         xiWtheta_ = ...
%     else
    xiWtheta = phiR/(2*sin(phiR))*dR;
%     end
end

%-------------- END CODE ---------------
end