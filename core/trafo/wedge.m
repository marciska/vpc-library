function y = wedge(x) %#codegen
%%WEDGE computes either u^ (x is 3-dim) or Vb^ (x is 6-dim)
%
% Syntax:
%
%   UWEDGE = WEDGE(U)
%       u is a vector of the form:
%               ⌈ a1 ⌉
%           u = | a2 |
%               ⌊ a3 ⌋
%       The wedge-transform of u is as follows:
%                ⌈  0  -a3  a2 ⌉
%           u^ = |  a3  0  -a1 |
%                ⌊ -a2  a1  0  ⌋
%
%   VBWEDGE = WEDGE(VB)
%       Vb is a vector of the form:
%                ⌈ v1 ⌉
%                | v2 |
%           Vb = | v3 |
%                | w1 |
%                | w2 |
%                ⌊ w3 ⌋
%       The wedge-transform of Vb is as follows:
%           Vb^ = ⌈   w^   v ⌉
%                 ⌊ 0 0 0  0 ⌋
%
%
% Inputs:
%
%   - U:  vector [size 3]
%   - VB: vector [size 6]
%
% Outputs:
%
%   - UWEDGE:  matrix [size 3x3]
%   - VBWEDGE: matrix [size 4x4]
%
% See also VEE.
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, March 2021
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

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
y = coder.nullcopy(zeros(ydim));
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
%%WEDGE computes u^
%
% Syntax:
%
%   UWEDGE = WEDGE(U)
%       u is a vector of the form:
%               ⌈ a1 ⌉
%           u = | a2 |
%               ⌊ a3 ⌋
%       The wedge-transform of u is as follows:
%                ⌈  0  -a3  a2 ⌉
%           u^ = |  a3  0  -a1 |
%                ⌊ -a2  a1  0  ⌋
%
%
% Inputs:
%
%   - U:  vector [size 3]
%
% Outputs:
%
%   - UWEDGE:  matrix [size 3x3]
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, March 2021
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

u_wedge = [  0  , -u(3),  u(2);
            u(3),   0  , -u(1);
           -u(2),  u(1),   0  ];

%-------------- END CODE ---------------
end

function Vb_wedge = vbwedge(Vb) %#codegen
%%VBWEDGE computes Vb^
%
% Syntax:
%
%   VBWEDGE = WEDGE(VB)
%       Vb is a vector of the form:
%                ⌈ v1 ⌉
%                | v2 |
%           Vb = | v3 |
%                | w1 |
%                | w2 |
%                ⌊ w3 ⌋
%       The wedge-transform of Vb is as follows:
%           Vb^ = ⌈   w^   v ⌉
%                 ⌊ 0 0 0  0 ⌋
%
%
% Inputs:
%
%   - VB: vector [size 6]
%
% Outputs:
%
%   - VBWEDGE: matrix [size 4x4]
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, March 2021
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

v = Vb(1:3);
omega = Vb(4:6);
Vb_wedge = [uwedge(omega), v(:);
             zeros(1,3),    0  ];

%-------------- END CODE ---------------
end