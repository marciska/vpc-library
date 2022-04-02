function sk_A = skew(A) %#codegen
%SKEW returns skew part of square matrix A
%
% Detailed Explanation:
%   none
%
% -----------
%
% Inputs:
%   - A: square matrix [size nxn]
%
% Outputs:
%   - sk_A: matrix (skew part of A) [size nxn]
%
% Example commands:
%   sk_A = skew([1 3; 2 10]);
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

sk_A = .5*(A-A');

%-------------- END CODE ---------------
end