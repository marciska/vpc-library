function y = wedge(x) %#codegen
%WEDGE computes either u^{wedge} (x is 3-dim) or Vb^{wedge} (x is 6-dim)
%
% Detailed Explanation:
%   * u is a vector of the form:
%         ⌈ a1 ⌉
%     u = | a2 |
%         ⌊ a3 ⌋
%     The wedge-transform of u is as follows:
%         ⌈  0  -a3  a2 ⌉
%     y = |  a3  0  -a1 |
%         ⌊ -a2  a1  0  ⌋
%   * Vb is a vector of the form:
%          ⌈ v1 ⌉
%          | v2 |
%     Vb = | v3 |
%          | w1 |
%          | w2 |
%          ⌊ w3 ⌋
%     The wedge-transform of Vb is as follows:
%     Vbwedge = ⌈ uwedge(w)  v ⌉
%               ⌊  0 0 0     0 ⌋
%
% -----------
%
% Inputs:
%   - x: vector [size 3 or size 6]
%
% Outputs:
%   - y: matrix [size 3x3 or size 4x4]
%
% Example commands:
%   uwedge  = wedge(ones(3,1));
%   vbwedge = wedge(ones(6,1));
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
% March 2021, Last modified: 2021-MAR-09
%
%------------- BEGIN CODE --------------

% Size of input vector (either 3 or 6)
n = length(x);

% Define size of output
ydim = 0;
switch n
    case 3
        ydim = 3;
    case 6
        ydim = 4;
    otherwise
        error('Input dimension is wrong!')
end

% Calculate wedge(x)
y = zeros(ydim);
switch n
    case 3
        y = uwedge(x);
    case 6
        y = vbwedge(x);
    otherwise
        error('Input dimension is wrong!')
end

%-------------- END CODE ---------------
end

function u_wedge = uwedge(u) %#codegen
%UWEDGE computes u^{wedge} for a 3-dimensional vector u
%
% Detailed Explanation:
%   u is a vector of the form:
%       ⌈ a1 ⌉
%   u = | a2 |
%       ⌊ a3 ⌋
%   The wedge-transform of u is as follows:
%             ⌈  0  -a3  a2 ⌉
%   u_wedge = |  a3  0  -a1 |
%             ⌊ -a2  a1  0  ⌋
%
% -----------
%
% Inputs:
%   - u: vector [size 3]
%
% Outputs:
%   - u_wedge: matrix [size 3x3]
%
% Example commands:
%   u_wedge = uwedge(ones(3,1));
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
% March 2021, Last modified: 2021-MAR-09
%
%------------- BEGIN CODE --------------

u_wedge = [  0  , -u(3),  u(2);
            u(3),   0  , -u(1);
           -u(2),  u(1),   0  ];

%-------------- END CODE ---------------
end

function vb_wedge = vbwedge(Vb) %#codegen
%VBWEDGE computes Vb^{wedge}
%
% Detailed Explanation:
%   Vb is a vector of the form:
%        ⌈ v1 ⌉
%        | v2 |
%   Vb = | v3 |
%        | w1 |
%        | w2 |
%        ⌊ w3 ⌋
%   The wedge-transform of Vb is as follows:
%   Vbwedge = ⌈ wedge(w)  v ⌉
%             ⌊  0 0 0    0 ⌋
%
% -----------
%
% Inputs:
%   - Vb: vector [size 6]
%
% Outputs:
%   - vb_wedge: matrix [size 4x4]
%
% Example commands:
%   vb_wedge = vbwedge(ones(6,1));
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
% March 2021, Last modified: 2021-MAR-09
%
%------------- BEGIN CODE --------------

v = Vb(1:3);
omega = Vb(4:6);
vb_wedge = [uwedge(omega), v(:);
             zeros(1,3),    0  ];

%-------------- END CODE ---------------
end