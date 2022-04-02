function f = vcamera(gco, lambda, po, camtype) %#codegen
%VCAMERA returns feature points in image coordinates based on pose g input.
%
% Detailed Explanation:
%   It is a virtual camera. That means, it can exactly calculate feature
%   point coordinates based on their (prior known) frame-position and their
%   current measured object pose.
%
% -----------
%
% Inputs:
%   - gco: pose matrix from object to camera [size 4x4]
%   - lambda: focal length [m] [scalar]
%   - po: feature point position matrix for n agents [m] [size nf x 3]:
%       po = [po1; po2; ...; ponf]
%   - camtype: camera ID representing a camera type [scalar]
%
% Outputs:
%   - f: feature points coordinates, vector [size (nf*2)]
%
% Example commands:
%   f = camera(eye(4));
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

% Put every po into homogeneous form [po; 1]
nf = size(po,1);
po_ = [po'; ones(1, nf)];

% Calculate visual measurement
f = coder.nullcopy(zeros(nf*2,1)); % Pre-allocate
switch camtype
    case 'Pinhole' % Pinhole camera
        % Transform feature points into camera frame Σc
        pc_ = gco*po_;

        % Form visual measurement vector f = [f1; f2; ...; fnf], with fi=[ui;vi]
        for i = 1:nf
            f(2*i-1:2*i) = lambda*pc_([1 3],i)/pc_(2,i);
        end
    otherwise % Unknown camera
        error('Unknown camera type!');
end

%-------------- END CODE ---------------
end