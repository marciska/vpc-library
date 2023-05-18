function kernel = SE3Hom(gv1, gv2, hyp, rho1, rho2) %#codegen
%%SE3 class of covariance kernel function for rigid motion.
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
%   - GV1: test pose input [size 4 x 4 x m ]
%   - GV2: test pose input [size 4 x 4 x m_]
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
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, 2023
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
[n1, n2, m] = size(gv1);
[~,  ~, m_] = size(gv2);
assert(n1==4 && n2 == 4,'Pose must be parsed in homogeneous form g = [R p; 0 1]!');

% extract hyperparameters
sigmafsq = hyp(1)^2;
lsq = hyp(2)^2;

% compute kernel
kernel = coder.nullcopy(zeros(m,m_));
for i = 1:m
    for j = 1:m_
        % split pose
        p1 = gv1(1:3,4,i);
        p2 = gv2(1:3,4,j);
        R1 = gv1(1:3,1:3,i);
        R2 = gv2(1:3,1:3,j);

        % compute distance
        dp = p1-p2;
%         distsq = rho1*dot(dp,dp) + rho2*(3 - (R1(:)')*R2(:))/2;
%         distsq = rho1*dot(dp,dp) + rho2*(3 - sum(sum(R1 .* R2)))/2;
        distsq = rho1*dot(dp,dp) + rho2*(norm(R1-R2,'fro')^2)/4;
        
        % compute kernel
        kernel(i,j) = sigmafsq * exp(-0.5*distsq/lsq);
    end
end

%-------------- END CODE ---------------
end
