function [mu, var] = gp_calc(xs, X, Y, kernel, hyp, sn, Kinv, A, varargin) %#codegen
%%GP_CALC returns the mean and variance of the GP Posterior distribution at
%%xs.
%
% Remarks:
%   For computational reasons, matrices Kinv and A can be parsed as
%   optional inputs. If they are not given, they will be comnputed
%   automatically.
%
% -----------
%
% Inputs:
%   - xs: test input [size 1xN]
%   - X: data set [size MxN]
%   - Y: data set [size MxQ]
%   - kernel: function-handler of kernel to be used [@functionHandler]
%   - sn: noise standard deviation on data Y [size Q]
%   - hyp: hyperparameters vector [size QxN]
% Optional Inputs:
%   - Kinv: Full kernel matrix [size MxMxQ]
%   - A: Precomputed data matrix, K^{-1}*Y [size MxQ]
%
% Outputs:
%   - mu: mean-vector of the multi-valued GP posterior [size 1xQ]
%   - var: variance-vector of the multi-valued GP posterior [size 1xQ]
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2023
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------
narginchk(5,inf)

% dimensions
Q = size(hyp,1);
M = size(Y,1);

% preset output dimension
mu = coder.nullcopy(zeros(1,Q));
var = coder.nullcopy(zeros(1,Q));

% check if data has been provided
if isempty(X) || isempty(Y)
    % default output is 0
    mu = zeros(1,Q);
    var = zeros(1,Q);
    return
end

% Calculate optional inputs if they have not been given
if nargin < 6 || isempty(sn)
    sn = zeros(Q,1);
end
if nargin < 7 || isempty(Kinv)
    Kinv = coder.nullcopy(zeros(M,M,Q));
    for q=1:Q
        L = chol(kernel(X,X,hyp(q,:),varargin{:}) + sn(q)^2*eye(M),"lower"); % A = L*L'
        Kinv(:,:,q) = (L\eye(M))'*(L\eye(M)); % A^-1 = L^-T*L^-1
        %Kinv(:,:,q) = inv(L)'*inv(L); % A^-1 = L^-T*L^-1
        %Kinv(:,:,q) = inv(kernel(X,X,hyp(q,:)) + sn(q)^2*eye(M));
    end
end
if nargin < 8 || isempty(A)
    A = coder.nullcopy(zeros(M,Q));
    for q=1:Q
        A(:,q) = Kinv(:,:,q) * Y(:,q);
    end
end

% Compute mean and variance
for q=1:Q
    % Kernel
    Kts = kernel(X,xs,hyp(q,:),varargin{:});
    Kst = Kts'; %kernel(xs,X,hyp(i,:),varargin{:});
    Kss = kernel(xs,xs,hyp(q,:),varargin{:});
    
    % mean
    mu(1,q) = Kst*A(:,q);
    var(1,q) = Kss - Kst*squeeze(Kinv(:,:,q))*Kts;
end

%-------------- END CODE ---------------
end