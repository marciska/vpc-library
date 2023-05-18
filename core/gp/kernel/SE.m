function kernel = SE(x1, x2, hyp) %#codegen
%SE represents the squared exponential covariance kernel class.
%
% Detailed Explanation:
%   none
%
% Remarks:
%   none
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
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2021
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% April 2021

%------------- BEGIN CODE --------------

% Dimensions
[m, ~] = size(x1);
[m_,~] = size(x2);

% Extract hyperparameters
sigmafsq = hyp(1)^2;
lsq = hyp(2)^2;

% Compute kernel
kernel = coder.nullcopy(zeros(m,m_));
for i = 1:m
    for j = 1:m_
        d = x1(i,:)-x2(j,:);
        r = dot(d,d)/lsq;
        kernel(i,j) = sigmafsq * exp(-r/2);
    end
end

%-------------- END CODE ---------------
end