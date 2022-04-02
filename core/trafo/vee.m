function y = vee(x) %#codegen
%VEE computes x^{v} (x can be 3x3 or 4x4)
%
% Detailed Explanation:
%   * x is a 3x3-matrix presented in wedge-transform:
%         ⌈  0  -a3  a2 ⌉
%     x = |  a3  0  -a1 |
%         ⌊ -a2  a1  0  ⌋
%     The vee-transform is the inverse operator of the wedge-transform, s.t.:
%         ⌈ a1 ⌉
%     y = | a2 |
%         ⌊ a3 ⌋
%   * x is a 4x4-matrix of the form:
%     x = ⌈ uwedge(w)  v ⌉
%         ⌊  0 0 0     0 ⌋
%     The vee-transform of x is as follows:
%         ⌈ v1 ⌉
%         | v2 |
%     y = | v3 |
%         | w1 |
%         | w2 |
%         ⌊ w3 ⌋
%
% -----------
%
% Inputs:
%   - x: matrix in wedge-transform [size 3x3 or 4x4]
%
% Outputs:
%   - y: vector [size 3 or 6]
%
% Example commands:
%   y = vee(...);
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
% January 2022
%
%------------- BEGIN CODE --------------

% size of input matrix (either 3 or 6)
n = size(x,1);

% calculate vee(x)
switch n
    case 3
        y = uvee(x);
    case 4
        y = Vbvee(x);
    otherwise
        error('Input dimension is wrong!')
end

%-------------- END CODE ---------------
end


function y = uvee(u) %#codegen
%VEE computes u^{v} (u is matrix of size 3x3)
%
% Detailed Explanation:
%   * u is a 3x3-matrix presented in wedge-transform:
%         ⌈  0  -a3  a2 ⌉
%     u = |  a3  0  -a1 |
%         ⌊ -a2  a1  0  ⌋
%     The vee-transform is the inverse operator of the wedge-transform, s.t.:
%         ⌈ a1 ⌉
%     y = | a2 |
%         ⌊ a3 ⌋
%
% -----------
%
% Inputs:
%   - u: matrix in wedge-transform [size 3x3]
%
% Outputs:
%   - y: vector [size 3]
%
% Example commands:
%   y = vee(...);
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
% January 2022
%
%------------- BEGIN CODE --------------

y = [u(3,2);
     u(1,3);
     u(2,1)];

%-------------- END CODE ---------------
end


function Vb = Vbvee(Vbwedge) %#codegen
%VEE computes Vbwedge^{v} (Vbwedge is matrix of size 4x4)
%
% Detailed Explanation:
%   * Vbwedge is a 4x4-matrix of the form:
%     Vbwedge = ⌈ uwedge(w)  v ⌉
%               ⌊  0 0 0     0 ⌋
%     The vee-transform of Vbwedge is as follows:
%          ⌈ v1 ⌉
%          | v2 |
%     Vb = | v3 |
%          | w1 |
%          | w2 |
%          ⌊ w3 ⌋
%
% -----------
%
% Inputs:
%   - Vbwedge: matrix in wedge-transform [size 4x4]
%
% Outputs:
%   - y: vector [size 6]
%
% Example commands:
%   y = vee(...);
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
% January 2022
%
%------------- BEGIN CODE --------------

v = Vbwedge(1:3,4);
w = uvee(Vbwedge(1:3,1:3));
Vb = [v; w];

%-------------- END CODE ---------------
end