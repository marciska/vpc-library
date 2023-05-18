function gv = check(g) %#codegen
%%CHECK returns vectorized form of homogeneous representation g.
%%In other words, g = [p xi*theta] is returned.
%
% Syntax:   GV = CHECK(G)
%
% Inputs:
%
%   - G: homogeneous representation / pose matrix [size 4x4]
%
% Outputs:
%
%   - GVEC: vector [size 6]
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, March 2021
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

if ndims(g) == 3
    m = size(g,3);
    gv = coder.nullcopy(zeros(m,6));
    % Check-transform all g
    for i=1:m
        gv(i,:) = check_(g(:,:,i));
    end
else
    % Check-transform g
    gv = check_(g);
end

%-------------- END CODE ---------------
end

function gv = check_(g) %#codegen
%%CHECK_ internal helper function for CHECK function.
%%Returns g = [p xi*theta].
%%See also CHECK.
%
% Syntax:   GV = CHECK_(G)
%
% Inputs:
%
%   - G: homogeneous representation / pose matrix [size 4x4]
%
% Outputs:
%
%   - GVEC: vector [size 6]
%
%
%   Editor: OMAINSKA Marco - Doctoral Student, Cybernetics
%               <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Property of: Fujita-Yamauchi Lab, The University of Tokyo, March 2021
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp

%------------- BEGIN CODE --------------

% Extract position & orientation from g
[R,p] = splitpose(g);

% Inverse-relation of exponential rotation matrix exp{xiwedge*theta}
% BUG: the following line has an issue when rot == pi. In practice, this
% will rarily create an issue because the numerics is never exactly pi in
% our simulation. Nevertheless, we should switch to 'rotm2axang' since it
% considers this corner case. Problem is, that rotm2axang creates some
% simulink errors about variable signals... we'll solve this at a later
% version of the VPC library
xitheta = vee(logR(R));
% axang = rotm2axang(R); % axang = [xi, theta]
% xitheta = axang(1:3).*axang(4);

% Form vector
gv = [p(:)' xitheta(:)'];

%-------------- END CODE ---------------
end