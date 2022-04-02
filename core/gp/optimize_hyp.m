function [hyp, sn] = optimize_hyp(X, Y, kernel, varargin) %#codegen
%OPTIMIZE_HYP returns the optimized Hyperparameters.
%
% Detailed Explanation:
%   Function commands as follows:
%       hyp = optimize_hyp(X, Y, kernel, ...)
%       [hyp, sn] = optimize_hyp(X, Y, kernel, ...)
%       _ = optimize_hyp(X, Y, kernel, standard_deviation)
%       _ = optimize_hyp(..., 'hyp', ones(data_length,xdim+1))
%       _ = optimize_hyp(..., 'maxsteps', 200)
%
% Remarks:
%   none
%
% -----------
%
% Inputs:
%   - X: data set [size MxN]
%   - Y: data set [size MxQ]
%   - kernel: function-handler of kernel to be used [@functionHandler]
%   - hyp: hyperparameters vector [size QxN]
%
% Outputs:
%   - mu: mean-vector of the multi-valued GP posterior [size Q]
%
% Example commands:
%   hyp = optimize_hyp([1 2 3; 4 5 6],[1 2; 7 8], @covSEard);
%   hyp = optimize_hyp([1 2 3; 4 5 6],[1 2; 7 8], @covSEard, 'maxsteps', 300);
%   hyp = optimize_hyp([1 2 3; 4 5 6],[1 2; 7 8], @covSEard, 1e-2*ones(2,1)); % adding noise
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

% add parsing parameters
p = inputParser;
addRequired(p,'X',@(x) ismatrix(x) && ~isempty(x));
addRequired(p,'Y',@(y) ismatrix(y) && ~isempty(y));
addRequired(p,'kernel',@(x) isa(x,'function_handle') || exist(x,'file'));

% Dimensions
[~,Q] = size(Y);
[~,N] = size(X);

% default parse parameters
default_sn = 1e-3*ones(Q,1);
if isa(kernel, 'string')
    kernel = str2func(kernel);
end
switch functions(kernel).file
    case functions(@covSE).file
        default_hyp = ones(Q,2);
    case functions(@covSEard).file
        default_hyp = ones(Q,N+1);
    otherwise
        error("Kernel not supported")
end

% add parsing parameters
addOptional(p,'sn',default_sn,@(x) isvector(x) && length(x)==Q)
addParameter(p,'hyp',default_hyp,@(x) ismatrix(x) && ~isempty(x));
addParameter(p,'maxsteps',100,@(x) isnumeric(x) && x > 0);

% parse inputs
parse(p,X,Y,kernel,varargin{:});
% disp(p.Results)

% optimize hyperparameters
hyp = p.Results.hyp;
sn = p.Results.sn;
for i=1:Q
    fprintf('## Optimize hyperparameters (%i/%i)\n',i,Q)
    hyper_opt = struct('cov',log(p.Results.hyp(i,:)),'lik',log(sn(i)));
    hyper_opt = minimize(hyper_opt, @gp, -p.Results.maxsteps, @infGaussLik, @meanZero, kernel, @likGauss, p.Results.X, p.Results.Y(:,i));
    
    % Copy optimized values to output in non-logarithmic space
    hyp(i,:) = exp(hyper_opt.cov);
    sn(i) = exp(hyper_opt.lik);
end


%-------------- END CODE ---------------
end