function R = lookRotation(a, b) %#codegen
%%LOOKROTATION computes rotation matrix R to rotate vector a onto vector b
%
% Detailed Explanation:
%   Rodrigues formula is used to calculate the rotation matrix. Please
%   refer to literature for more information.
%
% -----------
%
% Inputs:
%   - a: vector [size 3]
%   - b: vector [size 3]
%
% Outputs:
%   - R: matrix [size 3x3]
%
% Example commands:
%   R = lookRotation([1 0 0], [0 1 0]);
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
% May 2022

%------------- BEGIN CODE --------------

if all(a == 0) || all(b == 0)
    error('any of the parsed vectors is a zero vector');
end

% ensure a & b are unit vectors
a = a ./ norm(a);
b = b ./ norm(b);

if all(a == -b)
    error('parsed vectors look to opposite directions');
end

% important constants
v  = cross(a,b);
vw = wedge(v);   % skew symmetric cross product of v
%s2 = dot(v,v);   % sine-squared of angle
c  = dot(a,b);   % cosine of angle

% rodrigues formula
% R = eye(3) + vw + vw^2 .* (1-c)/s2; % can be simplified to below formula
R = eye(3) + vw + vw^2 .* 1/(1+c);

%------------- END CODE --------------
end