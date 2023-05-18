function kernel = SEARD(x1, x2, hyp) %#codegen
%SEARD squared-exponential (automatic relevance detection) kernel
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
%   - hyp: hyperparameters vector [size n+1]
%       - hyp[1]: signal ratio, sigmaf
%       - hyp[2:n+1]: lengthscales vector, l
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
% April 2021, Last modified: 2021-APR-26
%
%------------- BEGIN CODE --------------

% Dimensions
[m, n] = size(x1);
[m_,~] = size(x2);

% Pre-define kernel
kernel = coder.nullcopy(zeros(m,m_));

% Extract hyperparameters
sigmafsq = hyp(1)^2;
lsq = hyp(2:n+1).^2;
Lsq = diag(lsq);

% Compute SE-ARD kernel
for i = 1:m
    for j = 1:m_
        nu = (x1(i,:)-x2(j,:))';
        kernel(i,j) = sigmafsq * exp(-0.5*nu'*(Lsq\nu));
    end
end

%-------------- END CODE ---------------
end