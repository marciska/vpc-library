function [mu, cov] = gp_calc(x, X, Y, kernel, sn, hyp, Kinv, A) %#codegen
%GP_CALC returns the mean and covariance of the GP Posterior distribution at x
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
%   - x: test input [size N]
%   - X: data set [size MxN]
%   - Y: data set [size MxQ]
%   - kernel: function-handler of kernel to be used [@functionHandler]
%   - sigman: noise variance on data [size Q]
%   - hyp: hyperparameters vector [size QxN]
% Optional Inputs:
%   - Kinv: Full kernel matrix [size MxMxQ]
%   - A: Precomputed data matrix, K*Y [size MxQ]
%
% Outputs:
%   - mu: mean-vector of the multi-valued GP posterior [size Q]
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
%
%------------- BEGIN CODE --------------
narginchk(6,8)

% dimensions
Q = length(sn);
M = size(Y,1);

% preset output dimension
mu = coder.nullcopy(zeros(Q,1));
cov = coder.nullcopy(zeros(Q,1));

% bugfix: check if data has been provided
if isempty(X) || isempty(Y)
    % default output is 0
    mu = zeros(Q,1);
    cov = zeros(Q,1);
    return
end

% bugfix: make x always a row vector
x = x(:)';

% Calculate optional inputs if they have not been given
if nargin < 7
    Kinv = coder.nullcopy(zeros(M,M,Q));
    for i=1:Q
        Kinv(:,:,i) = inv(kernel(X,X,hyp(i,:)) + sn(i)^2*eye(M));
    end
end
if nargin < 8
    for i=1:Q
        A(:,i) = Kinv(:,:,i) * Y(:,i);
    end
    % A = pagemtimes(Kinv,Y); % not supported by MATLAB codegen
end

% Compute mean and cov
for i=1:Q
    % Kernel
    Kts = kernel(X,x,hyp(i,:));
    Kst = Kts'; %kernel(x,X,hyp(i,:));
    Kss = kernel(x,x,hyp(i,:));
    
    % mean
    mu(i) = Kst*A(:,i);
    cov(i) = Kss - Kst*squeeze(Kinv(:,:,i))*Kts;
end


%-------------- END CODE ---------------
end