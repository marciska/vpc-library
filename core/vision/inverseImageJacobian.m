function Jdagger = inverseImageJacobian(gco_bar, lambda, po) %#codegen
%INVERSEIMAGEJACOBIAN returns peusdo-inverse of image jacobian.
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - gco_bar: VMO-Estimate of pose matrix from object to camera [size 4x4]
%   - lambda: focal length [m] [scalar]
%   - po: feature point position matrix for n agents [m] [size nf x 3]:
%       po = [po1; po2; ...; ponf]
%
% Outputs:
%   - Jdagger: pseudo-inverse image jacobian [size 2nf x 2nf]
%
% Example commands:
%   Jdagger = imageJacobian(...);
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

% Extract pose
[Rco_bar,~] = splitpose(gco_bar);

% Put every po into homogeneous form [po; 1]
nf = size(po,1);
po_ = [po'; ones(1, nf)];

% Calculate estimated pc (feature points estimated by VMO)
pc_bar = gco_bar * po_;

% Calculate Image Jacobian
J = zeros(2*nf,6);
for i = 1:nf
    % Load i-th feature point (known & estimated)
    poi = po(i,:); % known
    pc_bar_i = pc_bar(:,i); % estimated
    
    % Calculate local image jacobian
%     Ji = [lambda/pc_bar_i(3),          0,         -lambda*pc_bar_i(1)/(pc_bar_i(3)^2); ...
%                   0,          lambda/pc_bar_i(3), -lambda*pc_bar_i(2)/(pc_bar_i(3)^2)];
    Ji = [lambda/pc_bar_i(2), -lambda*pc_bar_i(1)/(pc_bar_i(2)^2),         0; ...
                    0,        -lambda*pc_bar_i(3)/(pc_bar_i(2)^2), lambda/pc_bar_i(2)];
%     Ji = Ji * Rco_bar * [eye(3), -wedge(poi)];
    Ji = Ji * [Rco_bar, -Rco_bar*wedge(poi)];
    
    % Store to (global) image Jacobian
    J(2*i-1:2*i,:) = Ji;
end

% Transpose of image Jacobian
Jtrans = J';

% Pseudo-inverse of image Jacobian
Jdagger = (Jtrans*J)\Jtrans;

%-------------- END CODE ---------------
end