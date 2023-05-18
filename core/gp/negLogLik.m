function nloglik = negLogLik(X,Y,kernel,theta,logflag,varargin)
%%NEGLOGLIK computes the negative log-likelihood.
%
% Syntax:
%
%   NLOGLIK = NEGLOGLIK(X,Y,KERNEL,THETA)
%       Returns the negative log-likelihood at THETA for kernel KERNEL and
%       given data X, Y.
%
%   NLOGLIK = NEGLOGLIK(X,Y,KERNEL,LOGFLAG)
%       This setting determines if the parsed parameters at the
%       corresponding positions of the boolean vector LOGFLAG are in
%       logarithmic space. If so, the function will first transform them
%       back into normal space.
%           IMPORTANT: This setting can only be used for hyperparameters
%                      that can only be nonnegative.
%
%   NLOGLIK = NEGLOGLIK(X,Y,KERNEL,LOGFLAG,KERNELPARAMS)
%       The cell-array KERNELPARAMS will be parsed to the kernel KERNEL
%       in case the it has some parameters that are not classified as
%       hyperparameters.
%
%
% Inputs:
%
%   - X: input training data [size MxN]
%   - Y: output training data [size MxQ]
%   - KERNEL: kernel function [function handler or filename-string]
%   - HYP: GP kernel hyperparameter matrix [QxN_]
%   - SN: noise standard deviation vector/scalar [size Q / scalar]
%   - THETA: hyperparameter vector [size N_+1]
%       - THETA[1]: noise standard deviation SN [scalar]
%       - THETA[2:end]: kernel hyperparameters HYP [size N_]
%
% Optional Inputs:
%
%   - LOGFLAG: bool-vector determining if element of THETA is in log-space
%              [size N_+1]
%   - KERNELPARAMS: cell array of additional parameters that need to be
%                   parsed to the kernel
%
% Outputs:
%
%   - NLOGLIK: negative log-likelihood [scalar]
%
%
% See also OPTIMIZE_HYP.
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, May 2022
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------
algo = 1;
global encountered_non_positive_definite_kernel %#ok<GVMIS> 
narginchk(4,inf)
if nargin < 5 || isempty(logflag); logflag = false(size(theta)); end
assert(islogical(logflag),['islog MUST be of a logical datatype.' ...
    'Integer 0 or 1 will not work.']);

% transform inputs back if they are in log-space
logflag = boolean(logflag);
theta(logflag) = exp(theta(logflag));

% process inputs
sn = theta(1);
hyp = theta(2:end);
M = size(Y,1);

% gram matrix
K = kernel(X,X,hyp,varargin{:}) + sn^2*eye(M);

switch algo
    case 1 % cholesky
        try % try cholesky, else do fallback to normal matrix inversion
            % cholesky factorization of gram matrix: K = L*L'
            L = chol(K,'lower');

            % GP weight vector A
            A = L'\(L\Y);

            % negative log-likelihood
            nloglik = 0.5*dot(Y,A) + sum(log(diag(L))) + M/2*log(2*pi);
        catch ME
            % Note: sometimes the matrix is not perfectly positive definite
            % because of numerical precision.
            % This sometimes happens when the hyperparameters and noise are
            % badly set in one optimization step, thus rendering the
            % datapoints 'too close together'. When datapoints are too
            % close together, numerical inpreciseness results in the
            % mentioned error. 
            % This also means that the error can be prevented by setting
            % the hyperparameter/noise lower and upper bounds accordingly.
            % However, for now lets just do a fallback to normal matrix
            % inversion.
            encountered_non_positive_definite_kernel = true;

            % GP weight vector A
            A = K\Y;
            
            % negative log-likelihood
            nloglik = 0.5*dot(Y,A) + 0.5*det(K) + M/2*log(2*pi);
        end
    otherwise % normal matrix inversion
        % GP weight vector A
        A = K\Y;
        
        % negative log-likelihood
        nloglik = 0.5*dot(Y,A) + 0.5*det(K) + M/2*log(2*pi);
end

%-------------- END CODE ---------------
end
