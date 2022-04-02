function M = projectionSO3(M) %#codegen
%PROJECTIONSO3 returns the closest rotation matrix to matrix M in SO3.
%
% Detailed Explanation:
%   Projects an arbirary 3x3 matrix to the closest matrix in SO3.
%   Projection is given by formula:
%       Proj(M) = U*V'
%   with the terms from singular value decomposition
%       M = U*S*V'
%
% Remarks:
%   Due to integration, it might happen that Rot is not a real rotation
%   matrix anymore. But since Rot is "almost" a rotation matrix, it can be
%   projected into SO3 space.
%   If one does not do this, one might include orientation errors that get
%   bigger by integration over time. A possible way to check if a matrix is
%   not a real rotation anymore is to check the eigenvalues and
%   determinant.
%
% -----------
%
% Inputs:
%   - M: (arbitrary) rotation-like matrix [size 3x3]
%
% Outputs:
%   - M: rotation matrix [size 3x3]
%
% Example commands:
%   M = projection(eye(3));
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
% March 2021, Last modified: 2021-MAR-09
%
%------------- BEGIN CODE --------------

% Singular value decomposition
[U,~,V] = svd(M);

% Determinant of projected matrix (should be close to 1)
Rot_ = U*V.';
d = det(Rot_);

% Define identity matrix I, but set last element to d to make det(U*D*V.') exactly 1
I = diag([1 1 d]);

% Projection into SO3
M = U*I*V.';

%-------------- END CODE ---------------
end