function N = NAdj(gce) %#codegen
%NADJ calculates the matrix N
%
% Detailed Explanation:
%       ⌈            I_6                   0    ⌉
%   N = |                                       |
%       ⌊ -Ad(e^{-\xiwedge\theta_ce})     I_6   ⌋
%
% -----------
%
% Inputs:
%   - gce: pose control error matrix g [size 4x4]
%
% Outputs:
%   - N: Matrix [size 12x12]
%
% Example commands:
%   N = NAdj(eye(4))
%
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Review:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2020
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% Apr 2021, Last modified: 2021-APR-6
%
%------------- BEGIN CODE --------------

% Adjoint matrix
[Rce,~] = splitpose(gce);
Adj_ce = adjoint(Rce','rotation');

% N matrix

N = [  eye(6),  zeros(6);
      -Adj_ce,   eye(6) ];

%------------- END CODE --------------
end