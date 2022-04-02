function gvec = vec(g) %#codegen
%VEC returns vec-transform of pose g
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - g: pose matrix [size 4x4]
%
% Outputs:
%   - gvec: Vectorized transform of pose g [size 6]
%
% Example commands:
%   gvec = vec(mergepose(eye(3),[1 3 2]));
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

% Extract position & orientation from g
[R,p] = splitpose(g);

% Calculate skew(R)^v
skvee_R = vee(skew(R));

% Form vector
gvec = [p(:);
        skvee_R(:)];

%-------------- END CODE ---------------
end