function g_check = check(g) %#codegen
%CHECK returns check-transform of pose g
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - g: homogeneous representation / pose [size 4x4]
%
% Outputs:
%   - g_check: vector [size 6]
%
% Example commands:
%   g_check = check(eye(4));
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

if ndims(g) == 3
    m = size(g,3);
    g_check = coder.nullcopy(zeros(m,6));
    % Check-transform all g
    for i=1:m
        g_check(i,:) = check_(g(:,:,i));
    end
else
    % Check-transform g
    g_check = check_(g);
end

%-------------- END CODE ---------------
end

function g_check = check_(g) % internal check transform function

% Extract position & orientation from g
[R,p] = splitpose(g);

% Take inverse-relation of exponential rotation matrix exp{xiwedge*theta}
xitheta = vee(logR(R));

% Form vector
g_check = [p(:)' xitheta(:)'];

end