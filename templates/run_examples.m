%% run_examples.m
% *Summary:* Runs simple VMO+VPC experiments and plots the results
%
% -----------
%
% Editor:
%   OMAINSKA Marco - Doctoral Student, Cybernetics
%       <marcoomainska@g.ecc.u-tokyo.ac.jp>
% Supervisor:
%   YAMAUCHI Junya - Assistant Professor
%       <junya_yamauchi@ipc.i.u-tokyo.ac.jp>
%
% Property of: Fujita-Yamauchi Lab, University of Tokyo, 2021
% e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
% Website: https://www.scl.ipc.i.u-tokyo.ac.jp
% December 2021
%
% ------------- BEGIN CODE -------------

%% VMO

% simulate
simout_VMO = sim('VMO');

% 3D plot animation
figure('Name','VMO Animation','NumberTitle','off',...
    'Units','normalized','Position',[.05 .2 .4 .5]);
title('VMO Animation')
animate(simout_VMO);


%% VPC

% simulate
simout_VPC = sim('VPC');

% 3D plot animation
figure('Name','VPC Animation','NumberTitle','off',...
    'Units','normalized','Position',[.55 .2 .4 .5]);
title('VPC Animation')
animate(simout_VPC);


%% GP

% GP dataset
X = check(simout_VMO.gwo.signals.values);
% [~,X] = splitpose(simout_VMO.gwo.signals.values);
Y = simout_VMO.Vbwo.signals.values;

% bugfix: make constant signals time-varying signals
if length(simout_VMO.Vbwo) == 1
    Y = repmat(Y,size(X,1),1);
end

% reduce dataset
M = 20;
idx = randperm(size(X,1),M);
X = X(idx,:);
Y = Y(idx,:);
figure('Name','Dataset','NumberTitle','off',...
    'Units','normalized','Position',[.55 .2 .4 .5]);
plot3(X(:,1),X(:,2),X(:,3),'xr','MarkerSize',10,'LineWidth',2)

% add noise to dataset
sn = 1e-2;
Y = Y + sn^2.*gpml_randn(1,M,6);

% calculate GP hyperparameters
hyp = optimize_hyp(X,Y,@covSEard,sn*ones(6,1));
disp(hyp)


%% VMO + GP

% TODO: Set GP data in GP block

% simulate
simout_VMO_GP = sim('VMO_GP');

% 3D plot animation
figure('Name','VMO+GP Animation','NumberTitle','off',...
    'Units','normalized','Position',[.55 .2 .4 .5]);
if exist('X','var')
    plot3(X(:,1),X(:,2),X(:,3),'xr','MarkerSize',10,'LineWidth',2)
end
animate(simout_VMO_GP);


%% VPC + GP

% TODO: Set GP data in GP block

% simulate
simout_VPC_GP = sim('VPC_GP');

% 3D plot animation
figure('Name','VPC+GP Animation','NumberTitle','off',...
    'Units','normalized','Position',[.55 .2 .4 .5]);
if exist('X','var')
    plot3(X(:,1),X(:,2),X(:,3),'xr','MarkerSize',10,'LineWidth',2)
end
animate(simout_VPC_GP);


%% VMO + onlineGP: time-triggered

% % simulate
% simout_VMO_GP_triggerTime = sim('VMO_onlineGP_triggerTime');
% 
% % 3D plot animation
% figure('Name','VMO+onlineGP (Time-Trigger) Animation','NumberTitle','off',...
%     'Units','normalized','Position',[.55 .2 .4 .5]);
% animate(simout_VMO_GP_triggerTime);


%% VPC + onlineGP: time-triggered

% % simulate
% simout_VPC_GP_triggerTime = sim('VPC_onlineGP_triggerTime');
% 
% % 3D plot animation
% figure('Name','VPC+onlineGP (Time-Trigger) Animation','NumberTitle','off',...
%     'Units','normalized','Position',[.55 .2 .4 .5]);
% animate(simout_VPC_GP_triggerTime);


%% VMO + onlineGP: Lyapunov-triggered

% % simulate
% simout_VMO_GP_triggerLyapunov = sim('VMO_onlineGP_triggerLyapunov');
% 
% % 3D plot animation
% figure('Name','VMO+onlineGP (Lyapunov-Trigger) Animation','NumberTitle','off',...
%     'Units','normalized','Position',[.55 .2 .4 .5]);
% animate(simout_VMO_GP_triggerLyapunov);


%% VPC + onlineGP: Lyapunov-triggered

% % simulate
% simout_VPC_GP_triggerLyapunov = sim('VPC_onlineGP_triggerLyapunov');
% 
% % 3D plot animation
% figure('Name','VPC+onlineGP (Lyapunov-Trigger) Animation','NumberTitle','off',...
%     'Units','normalized','Position',[.55 .2 .4 .5]);
% animate(simout_VPC_GP_triggerLyapunov);
