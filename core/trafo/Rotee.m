function Ree = Rotee(x) %#codegen
%ROTEE calculates the rotation estimation error Ree from the estimation
%error ee or directlty from skRv = sk(e^{xi^{∨}*θ_ee})^{∨}
%
% Detailed Explanation:
%   Under the assumption that the estimation error angle is bounded by
%   |θ_ee| ≤ π/2, the estimation rotation error R = e^{xi^{∨}*θ_ee}
%   can be calculated by its vectorized skew-form skRv = sk(R)^{∨} as:
%       R = asin(||skRv||)/||skRv|| * skRv
%
% -----------
%
% Inputs:
%   - x: Either skRv or ee can be passed
%       - ee: estimation error vector [size 6]
%       - skRv: vee-transformed skew rotation error matrix [size 3]
%
% Outputs:
%   - Ree: matrix [size 3x3]
%
% Example commands:
%   Ree = Rotee(zeros(6,1));
%   Ree = Rotee(zeros(3,1));
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita Lab, University of Tokyo, 2021
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% December 2021
%
%------------- BEGIN CODE --------------

% parse inputs
switch length(x)
    case 3
        skRv = x;
    case 6
        skRv = x(4:6);
    otherwise
        error('Parsed input is neither ee nor skRv!')
end

% Compute rotation estimation error
if nnz(skRv) == 0 % if input is 0-vector
    Ree = eye(3); % default rotation matrix (no rotation estimation error)
else
    % Calculate ||sk(e^{\xi\theta_{ee}})^{v}||
    abs_skRv = min(norm(skRv),1);
    
    % Compute axis-angle pair
    xi_theta_ee = asin(abs_skRv)/abs_skRv * skRv;

    % Compute rotation estimation error matrix
    Ree = expm(wedge(xi_theta_ee));
end

%-------------- END CODE ---------------
end