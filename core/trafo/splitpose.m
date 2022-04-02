function [R,p] = splitpose(g) %#codegen
%SPLITPOSE returns Rotation R and position p from pose g
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - g: Homogeneous representation / Pose matrix [size 4x4]
%
% Outputs:
%   - R: Rotation matrix [size 3x3]
%   - p: Position vector [size 3]
%
% Example commands:
%   [R,p] = splitpose(eye(4));
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
% March 2021, Last modified: 2021-MAR-08
%
%------------- BEGIN CODE --------------

dim = ndims(g);
switch dim
    case 2
        R = g(1:3,1:3,:);
        p = g(1:3,4,:);
    case 3
        M = size(g,3);
        R = zeros(3,3,M);
        p = zeros(M,3);
        for i = 1:M
            R(:,:,i) = g(1:3,1:3,i);
            p(i,:) = g(1:3,4,i)';
        end
    otherwise
        error("Did not provide pose data. Couldn't split data.")
end

%------------- END CODE --------------
end