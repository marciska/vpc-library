function g = unvec(x) %#codegen
%UNVEC is the back-transform of the vec-transform vec(g)
%
% Note:
%   This back-transformation does only hold under the assumption that the
%   angle is bounded as |θ| ≤ π/2. Only under this assumption the following
%   formula holds:
%       R = asin(||skRv||)/||skRv|| * skRv
%
% -----------
%
% Inputs:
%   - x: Vectorized transform of pose g [size 6]
%
% Outputs:
%   - g: pose matrix g [size 4x4]
%
% Example commands:
%   g = unvec([1 -2 1 0 0 1]);
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
% April 2022
%
%------------- BEGIN CODE --------------

% Extract position & orientation from x
p = x(1:3);
skRv = x(4:6);

% Compute rotation matrix
R = eye(3); % default rotation matrix (no rotation)
if nnz(skRv) ~= 0 % if input is NOT a 0-vector
    % Calculate ||sk(e^{\xi\theta_{ee}})^{v}||
    abs_skRv = min(norm(skRv),1);
    % Compute axis-angle pair
    xi_theta = asin(abs_skRv)/abs_skRv * skRv;
    % Compute rotation matrix
    R = expm(wedge(xi_theta));
end

% Pose g
g = mergepose(R,p);

%-------------- END CODE ---------------
end