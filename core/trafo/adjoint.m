function Ad = adjoint(data, adjointtype) %#codegen
%ADJOINT calculates adjoint matrix
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - data: data adjoint matrix is calculated for. Can be one of these:
%       * pose matrix g [size 4x4]
%       * position vector p [size 3]
%       * orientation matrix R [size 3x3]
%   - adjointtype: type adjoint matrix is calculated for. Can be one of:
%       * 'pose'
%       * 'position'
%       * 'orientation'
%
% Outputs:
%   - Ad: Adjoint matrix [size 6x6]
%
% Example commands:
%   Ad = adjoint(eye(4),'pose')
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

% Assign position & orientation based on data provided
p = zeros(3,1);
R = zeros(3);
switch adjointtype
    case 'pose'
        [R,p] = splitpose(data);
    case 'position'
        p = data;
    case {'orientation','rotation'}
        R = data;
    otherwise
        error('Adjoint type not supported');
end

% Adjoint matrix
Ad = [   R,     wedge(p)*R;
      zeros(3),      R    ];

%------------- END CODE --------------
end