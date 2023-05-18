function kernel = SE3Axang(gv1, gv2, hyp, rho1, rho2) %#codegen
%%SE3 class of covariance kernel function for rigid motion.
%
% Kernel has been proposed in the following paper:
%   Lang, M. & Hirche, S. (2017). Computationally Efficient Rigid-Body
%   Gaussian Process for Motion Dynamics. IEEE Robotics and Automation
%   Letters, 2(3), 1601–1608. https://doi.org/10.1109/LRA.2017.2677469
%
% Syntax:
%
%   KERNEL = SE3(GV1, GV2, HYP)
%       Computes the SE(3) kernel for poses GV1 and GV2 based on the kernel
%       hyperparameters HYP and [ρ1 ρ2]=[1/√2 1/√2].
%
%   KERNEL = SE3(GV1, GV2, HYP, RHO1, RHO2)
%       Computes the SE(3) kernel for poses GV1 and GV2 based on the kernel
%       hyperparameters HYP and [ρ1 ρ2]=[RHO1 RHO2]. Note that vector
%       [ρ1 ρ2] will be normalized beforehand.
%
%
% Inputs:
%
%   - GV1: test pose input [size m  x 6]
%   - GV2: test pose input [size m_ x 6]
%   - HYP: hyperparameters vector [size 2]
%       - hyp[1]: signal ratio σ
%       - hyp[2]: lengthscale l
%   - RHO1: weight for rotation distance [scalar]
%   - RHO2: weight for translation distance [scalar]
%
% Outputs:
%
%   - KERNEL: kernel matrix [size mxm_]
%
% See also SE, SEARD, MATERN, MATERNARD.
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, May 2022
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

% pre-process inputs
narginchk(3,5)
if nargin < 4
    rho1 = 1/2;
    rho2 = 1/2;
elseif nargin < 5
    rho2 = 1 - rho1;
end
rho_norm = rho1 + rho2;
rho1 = rho1/rho_norm; rho2 = rho2/rho_norm;

% dimensions
[m, n] = size(gv1);
[m_,~] = size(gv2);
assert(n==6,'Pose must be parsed in vector form g = [p xitheta]!');

% extract hyperparameters
sigmafsq = hyp(1)^2;
lsq = hyp(2)^2;

% compute kernel
kernel = coder.nullcopy(zeros(m,m_));
for i = 1:m
    for j = 1:m_
        % split pose
        p1 = gv1(i,1:3); xitheta1 = gv1(i,4:6);
        p2 = gv2(j,1:3); xitheta2 = gv2(j,4:6);

        % compute distance
        dp = p1-p2;
        distsq = rho1*dot(dp,dp) + rho2*darc(xitheta1,xitheta2)^2;
        
        % compute kernel
        kernel(i,j) = sigmafsq * exp(-0.5*distsq/lsq);
    end
end

%-------------- END CODE ---------------
end

function d = darc(xitheta1,xitheta2)
%%DARC distance function metric for axis-angle rotations 
%
% Syntax:   D = DARC(XITHETA1,XITHETA2)
%
% Inputs:
%   - XITHETA1: axis-angle rotation [size 3]
%   - XITHETA2: axis-angle rotation [size 3]
%
% Outputs:
%   - D: distance [scalar]

%------------- BEGIN CODE --------------

% pre-process inputs
[xi1, theta1] = splitAxisAngle(xitheta1);
[xi2, theta2] = splitAxisAngle(xitheta2);

% calculate distance
% Note: Sometimes x is larger than 1, which is a problem for acos().
% This error comes just from numerics, so we saturate x as: -1 <= x <= 1
x = abs( cos(theta1/2)*cos(theta2/2) ...
        + sin(theta1/2)*sin(theta2/2)*dot(xi1,xi2) );
x = max(-1, min(x,1));
d = 2*acos(x);

%-------------- END CODE ---------------
end

function [xi, theta] = splitAxisAngle(xitheta)
%%SPLITAXISANGLE splits the axis-angle composed vector in axis and angle
%
% Syntax:   [XI,THETA] = SPLITAXISANGLE(XITHETA)
%
% Inputs:
%   - XITHETA: axis-angle rotation [size 3]
%
% Outputs:
%   - XI: rotation axis [size 3]
%   - THETA: rotation angle [scalar]

%------------- BEGIN CODE --------------

theta = norm(xitheta,2);
if theta == 0 % corner case
    xi = [0 0 1]; % default axis
else
    xi = xitheta./theta;
end

%-------------- END CODE ---------------
end