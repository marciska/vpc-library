function [hyp, sn, nloglik] = optimize_hyp(X,Y,kernel,varargin)
%%OPTIMIZE_HYP returns the optimized hyperparameters of a Gaussian
%%Process by minimizing the negative log-likelihood.
%
% Syntax:
%
%   [HYP,SN,NLOGLIK] = OPTIMIZE_HYP(X,Y,KERNEL,[HYP])
%       Optimizes the hyperparameters HYP of KERNEL for the given data X
%       and Y. If KERNEL is a function from this library, the initial HYP
%       matrix does not need to be given. The optimized noise standard
%       deviation SN and negative log-likelihood NLOGLIK are also returned.
%
%   [___] = OPTIMIZE_HYP(___,'sn0',SN)
%       Starts the optimization from initial standard deviation SN.
%
%   [___] = OPTIMIZE_HYP(___,'kernelparams',KERNELPARAMS)
%       The cell-array KERNELPARAMS will be parsed to the KERNEL function
%       in case the kernel has some parameters that are not classified as
%       hyperparameters to be optimized.
%
%   [___] = OPTIMIZE_HYP(___,'fixedHyper',FHYP)
%       If some hyperparameters shall be kept constant during optimization,
%       parse in this boolean matrix FHYP of parameter positions.
%
%   [___] = OPTIMIZE_HYP(___,'fixedNoise',FSNFLAG)
%       The flag FSNFLAG determines if the standard deviation SN should be
%       fixed during optimization.
%
%   [___] = OPTIMIZE_HYP(___,'lowerBound',LB)
%       LB is a vector of lower-bounds of the hyperparameters that
%       optimization has to fulfill.
%
%   [___] = OPTIMIZE_HYP(___,'upperBound',UB)
%       UB is a vector of upper-bounds of the hyperparameters that
%       optimization has to fulfill.
%
%   [___] = OPTIMIZE_HYP(___,'solveInLogSpace',LOGFLAG)
%       This setting will optimize the given parameter-positions LOGFLAG in
%       logarithmic space. Usually, the optimization in log-space is faster
%       and sometimes comes to a better solution, but experimenting with
%       this setting on a case-by-case is recommended.
%           IMPORTANT: This setting can only be used for hyperparameters
%                      that can only be nonnegative.
%
%   [___] = OPTIMIZE_HYP(___,'optimoptions',OPTIONS)
%       Cell- or struct-array of additional key-value pairs that should be
%       parsed to the optimizer-function FMINCON.
%       See also OPTIMOPTIONS, FMINCON.
%
%   [___] = OPTIMIZE_HYP(___,'verbose',VERBOSE)
%       VERBOSE flag indicates if all optimization-iterations and technical
%       details shall be printed to the MATLAB console.
%
%
% Inputs:
% 
%   - X: input training data [size MxN]
%   - Y: output training data [size MxQ]
%   - KERNEL: kernel function [function handler or filename-string]
%   - HYP: GP kernel hyperparameter matrix [QxN_]
%   - SN: noise standard deviation vector/scalar [size Q / scalar]
%
% Optional Inputs [key-value pair]:
%
%   - KERNELPARAMS: cell array of additional parameters that need to be
%                   parsed to the kernel
%   - FHYP: bool-matrix determining which hyperparameters shall be kept
%           constant during optimization [size QxN_]
%   - FSNFLAG: flag determining if the noise level shall be kept constant
%              during optimization [size Q]
%   - LB: vector of lower-bounds for hyperparameters [size N_]
%   - UB: vector of upper-bounds for hyperparameters [size N_]
%   - LOGFLAG: bool-vector determining which hyperparameters shall be
%              optimized in logarithmic space [bool, size N_]
%   - OPTIONS: cell or struct array of parameters that will be parsed to
%              the fmincon-function (the optimizer)
%   - VERBOSE: flag indicating if iterations shall be printed [bool]
%
% Outputs:
%
%   - hyp: optimized GP kernel hyperparameters [size QxN_]
%   - sn: optimized noise standard deviation vector [size Q]
%   - nloglik: optimized negative log-likelihood at solution [size Q]
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, May 2022
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

%% Parsing inputs

% pre-process some important inputs
% validateattributes(X,{'double'},{'nonempty','2d','finite','nonnan'},'','X')
validateattributes(Y,{'double'},{'nonempty','2d','finite','nonnan'},'','Y')
assert(isa(kernel,'function_handle') || exist(kernel,'file'),...
  'Invalid kernel given. Must be a function handle or string of filename');
if isa(kernel, 'string'); kernel = str2func(kernel); end
assert(abs(nargin(kernel))>=3,['It seems the parsed kernel function is' ...
    ' invalid. Check if it follows this pattern: @kernel(X,Y,hyp, ...)'])
D = ndims(X);
if D==2
    N = size(X,2);
end
Q = size(Y,2);

% default values for optional parameters
default_hyp = [];
default_solveInLogSpace = false;
default_lb = [];
default_ub = [];
switch functions(kernel).file
    case {functions(@SE).file, functions(@Matern).file, ...
            functions(@SE3Axang).file, functions(@SE3Hom).file}
        default_hyp = ones(Q,2);
        default_lb = 1e-3.*ones(Q,2);
        default_ub = 1e3.*ones(Q,2);
        default_solveInLogSpace = true;
    case {functions(@SEARD).file, functions(@MaternARD).file}
        default_hyp = ones(Q,N+1);
        default_lb = 1e-3.*ones(Q,N+1);
        default_ub = 1e3.*ones(Q,N+1);
        default_solveInLogSpace = true;
end
default_sn = 1e-1;
lb_sn = 1e-3;
ub_sn = 1e1;

% parse inputs
p = inputParser;
switch functions(kernel).file
    case {functions(@SE).file, functions(@SEARD).file, ...
          functions(@Matern).file, functions(@MaternARD).file, ...
          functions(@SE3Axang).file, functions(@SE3Hom).file}
        addOptional(p,'hyp',default_hyp,@(x)validateattributes(x,...
            {'numeric'},{'2d','finite','nonnan','positive'}));
    otherwise
        assert(~isstring(varargin{1}) && ~ischar(varargin{1}), ...
            ['Cannot infer length of hyperparameter vector for given ' ...
            'kernel. It must be specified.']);
        addRequired(p,'hyp',@(x)validateattributes(x,{'numeric'},...
            {'2d','finite','nonnan'}));
end
addParameter(p,'sn0',default_sn,@(x)validateattributes(x,...
    {'numeric'},{'finite','>=',lb_sn,'<=',ub_sn}));
addParameter(p,'kernelparams',{},@(x)validateattributes(x,{'cell'},{}));
addParameter(p,'fixedHyper',false,@(x)validateattributes(x,...
            {'logical','numeric'},{'2d','nonempty','binary'}));
addParameter(p,'fixedNoise',false,@(x)validateattributes(x,...
            {'logical','numeric'},{'scalar','binary'}));
addParameter(p,'lowerBound',default_lb,@(x)validateattributes(x,...
    {'numeric'},{'nonnan'}));
addParameter(p,'upperBound',default_ub,@(x)validateattributes(x,...
    {'numeric'},{'nonnan'}));
addParameter(p,'solveInLogSpace',default_solveInLogSpace,...
    @(x)validateattributes(x,{'logical','numeric'},{'vector','binary'}));
addParameter(p,'optimoptions',{},@(x)validateattributes(x,...
    {'cell','struct'},{}));
addParameter(p,'verbose',true,@(x)validateattributes(x,...
    {'logical','numeric'},{'scalar','binary'}));
parse(p,varargin{:});

% process parsed parameters
hyp = p.Results.hyp;
if isempty(hyp); hyp = default_hyp; end
assert(~isempty(hyp),['A kernel was parsed without prior info about ' ...
    'the size of the hyperparameter vector. Please specify an initial ' ...
    'hyperparameter vector to be optimized'])
fixedHyper = logical(p.Results.fixedHyper);
if length(fixedHyper)==1
    fixedHyper = logical(fixedHyper.*ones(Q,size(hyp,2)));
elseif all([1 size(hyp,2)]==size(fixedHyper))
    fixedHyper = repmat(fixedHyper,Q,1);
elseif all([Q 1]==size(fixedHyper))
    fixedHyper = repmat(fixedHyper,1,size(hyp,2));
end
assert(all(size(fixedHyper)==size(hyp)),['Error in dimension length of' ...
    'fixedHyper: %s~=%s'],mat2str(size(fixedHyper)),mat2str(size(hyp)))
sn = p.Results.sn0; sn = sn(:);
if isempty(sn); sn = default_sn; end
if length(sn)==1; sn = sn.*ones(Q,1); end
assert(length(sn)==Q,'Expected sn to be an array of %i elements',Q);
fixedNoise = logical(p.Results.fixedNoise(:));
% if length(fixedNoise)==1; fixedNoise = logical(fixedNoise.*ones(Q,1)); end
% assert(length(fixedNoise)==Q,['Expected fixedNoise to be an array of ' ...
%     '%i elements'],Q)
solveInLogSpace = logical(p.Results.solveInLogSpace(:)');
if length(solveInLogSpace)==1; ...
       solveInLogSpace = logical(solveInLogSpace.*ones(1,size(hyp,2))); end
assert(length(solveInLogSpace)==size(hyp,2),['Expected solveInLogSpace' ...
    ' to be an array of %i elements'],size(hyp,2))
solveInLogSpace = [true solveInLogSpace]; % sn is always optimized in log
lb_hyp = p.Results.lowerBound;
if all([1 size(hyp,2)]==size(lb_hyp)); lb_hyp = repmat(lb_hyp,Q,1); end
if ~isempty(lb_hyp)
    assert(all(size(lb_hyp)==size(hyp)),['Expected lowerBound to be an' ...
        ' array of %i elements'],size(hyp,2))
    assert(all(lb_hyp(:,solveInLogSpace(2:end))>=0,'all'),['Some ' ...
        'hyperparameters were requested to be optimized in log-space, ' ...
        'but their lower bound is negative, which is not possible. ' ...
        'Either remove the solveInLogSpace flag vor the corresponding ' ...
        'hyperparameters, or set their lower bound to be nonnegative.'])
end
ub_hyp = p.Results.upperBound;
if all([1 size(hyp,2)]==size(ub_hyp)); ub_hyp = repmat(ub_hyp,Q,1); end
if ~isempty(ub_hyp)
    assert(all(size(ub_hyp)==size(hyp)),['Expected upperBound to be an ' ...
        'array of %i elements'],size(hyp,2))
    assert(all(ub_hyp(:,solveInLogSpace(2:end))>=0,'all'),['Some ' ...
        'hyperparameters were requested to be optimized in log-space, ' ...
        'but their lower bound is negative, which is not possible. ' ...
        'Either remove the solveInLogSpace flag vor the corresponding ' ...
        'hyperparameters, or set their lower bound to be nonnegative.'])
end
lb_sn = lb_sn.*ones(Q,1);
ub_sn = ub_sn.*ones(Q,1);
kernelparams = p.Results.kernelparams;
options = p.Results.optimoptions;
if isstruct(options); options = namedargs2cell(options); end
verbose = p.Results.verbose;


%% Optimization
global encountered_non_positive_definite_kernel %#ok<GVMIS>
encountered_non_positive_definite_kernel = false;

% options
if verbose; displaytype = 'iter'; else; displaytype = 'none'; end
options = optimoptions(@fmincon,'Display',displaytype,options{:});

% loop over all GP outputs
if verbose; fprintf('Starting hyperparameter optimization.\n'); tic; end
nloglik = coder.nullcopy(zeros(Q,1));
for q = 1:Q
    % initial starting point
    theta0 = [sn(q), hyp(q,:)];
    n = length(theta0);
    
    % minimization problem
    fun = @(x)negLogLik(X,Y(:,q),kernel,x,solveInLogSpace,kernelparams{:});

    % constraints
    fixedTheta = [fixedNoise fixedHyper(q,:)];
    Aeq = []; beq = [];
    if any(fixedTheta)
        Aeq = zeros(n,n); beq = zeros(n,1);
        Aeq(fixedTheta,fixedTheta) = 1;
        beq(fixedTheta) = theta0(fixedTheta);
    end
    lb = -inf(1,n); ub = inf(1,n);
    if ~isempty(lb_sn); lb(1) = lb_sn(q); end
    if ~isempty(lb_hyp); lb(2:n) = lb_hyp(q,:); end
    if all(lb==-inf); lb = []; end
    if ~isempty(ub_sn); ub(1) = ub_sn(q); end
    if ~isempty(ub_hyp); ub(2:n) = ub_hyp(q,:); end
    if all(ub==inf); ub = []; end
    
    % transform constraints and IC to log-space if requested
    theta0(solveInLogSpace) = log(theta0(solveInLogSpace));
    if ~isempty(beq)
        logAndFixed = solveInLogSpace & fixedTheta;
        beq(logAndFixed) = log(beq(logAndFixed)); %#ok<AGROW>
                                % some linter error I guess. Because beq
                                % does not change size on each iteration.
    end
    if ~isempty(lb); lb(solveInLogSpace) = log(lb(solveInLogSpace)); end
    if ~isempty(ub); ub(solveInLogSpace) = log(ub(solveInLogSpace)); end
    
    % minimize
    if verbose; fprintf('# (%i/%i) GP optimization\n',q,Q); end
    [theta, nloglik_] = fmincon(fun,theta0,[],[],Aeq,beq,lb,ub,[],options);
    theta(solveInLogSpace) = exp(theta(solveInLogSpace));
    if verbose; fprintf('Done with nloglik=%g\n\n',nloglik_); end

    % set outputs
    nloglik(q) = nloglik_;
    sn(q) = theta(1);
    hyp(q,:) = theta(2:end);
end

% print optimized hyperparameters
if verbose
    fprintf('~~~~ Optimization terminated after t=%.2f sec.\n',toc);
    disp('Hyperparameters [optimized]:')
    disp(hyp)
    disp('Standard Deviation [optimized]:')
    disp(sn(:)')
    disp('Negative Log-Likelihood [optimized]:')
    disp(nloglik(:)')
end
if encountered_non_positive_definite_kernel
    warning(['Encountered a non-positive-definite matrix during ' ...
        'optimization. You might want to set the hyperparameter and ' ...
        'noise lower/upper-bounds so that this error does not occur. ' ...
        'Usually, this error results in degraded GP performance and ' ...
        'should be avoided.']); %#ok<UNRCH> 
end
clear encountered_non_positive_definite_kernel

%-------------- END CODE ---------------
end
