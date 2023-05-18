function kernel = Matern(x1, x2, hyp, p) %#codegen
%%Matérn class of covariance function.
%
% Detailed Explanation:
%   Specialized mattern covariance kernel class for nu = p + 1/2. The
%   function accepts only p ∈ {0, 1, 2} which result in the following
%   simplified representations of the mattern kernel:
%       p=0: k(r) = σ * exp(-r)
%       p=1: k(r) = σ * (1 + √(3)r) * exp(-√(3)r)
%       p=2: k(r) = σ * (1 + √(5)r + (5/3)r²) * exp(-√(5)r)
%   The signal ratio is represented by σ and the distance is defined as
%   r = |x1-x2|/l² with l the lengthscale. Both σ and l are
%   hyperparameters. p can be given by the function call, whereas if it was
%   not given then the default value is p = 2.
%
% -----------
%
% Inputs:
%   - x1: test input [size m  x n]
%   - x2: test input [size m_ x n]
%   - hyp: hyperparameters vector [size 2]
%       - hyp[1]: signal ratio σ
%       - hyp[2]: lengthscale l
%
% Outputs:
%   - kernel: kernel matrix [size mxm_]
%
% Example commands:
%   none
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2022
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% May 2022

%------------- BEGIN CODE --------------
narginchk(3,4)
if nargin < 4
    p = 2;
end

% Dimensions
[m, ~] = size(x1);
[m_,~] = size(x2);

% Extract hyperparameters
sigmafsq = hyp(1)^2;
l = hyp(2);

% Compute kernel
kernel = coder.nullcopy(zeros(m,m_));
for i = 1:m
    for j = 1:m_
        r = norm(x1(i,:)-x2(j,:),2)/l;
        switch p
            case 0
                kernel(i,j) = sigmafsq * exp(-r);
            case 1
                kernel(i,j) = sigmafsq * (1+sqrt(3)*r) * exp(-sqrt(3)*r);
            case 2
                kernel(i,j) = sigmafsq * (1+sqrt(5)*r + (5/3)*r^2) * exp(-sqrt(5)*r);
            otherwise
                error('Invalid p given (p=%g). Valid values are {0,1,2}.',p)
        end
    end
end

%-------------- END CODE ---------------
end