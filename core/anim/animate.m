function animate(varargin)
%ANIMATE VPC simulation.
%
% Detailed Explanation:
%   TODO
%
% Function parameter calls:
%   TODO
%
% -----------
%
% Inputs:
%   - ax: axes-object
%   - simout: simulation-output from SIMULINK
%   - TODO
%
% Example commands:
%   - animate(simout)
%   - animate(simout,'bound',[-4,4;-4,4;-4,4])
%   - animate(simout,'fps',60)
%   - animate(simout,'fps',60,'recorder',struct('file','out.mp4','profile','MPEG-4','quality',100))
%   - animate(...
%       {out.gwo.signals.values,agentdesign('target','show_trajectory','on')},...
%       {out.gwc.signals.values,agentdesign('agent'),out.gcobar.signals.values,agentdesign('estimate')})
%   - animate(...
%       {out.gwo.signals.values,agentdesign('target','show_trajectory','on')},...
%       {out.gwc.signals.values,agentdesign('agent'),out.gcobar.signals.values,agentdesign('estimate')},...
%       'communication',[1 2])
%   - animate(...
%       out.tout,{out.gwo.signals.values,agentdesign('target','show_trajectory','on')},...
%       'fps',60,'recorder',struct('file','out.mp4','profile','MPEG-4','quality',100))
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
% December 2021
%
% ------------- BEGIN CODE --------------

% TODO
% - axis-relation: [1 1 1]  --> even if 1 axis is differently bounded, it
%       should not distort the graphics
% - make gradient-trajectory position-dependent, not time-index-dependent


%% Process inputs

% Input Channel iterator
k = 1;

% Check if axes-object has been provided
if isa(varargin{k}, 'matlab.graphics.axis.Axes')
    ax = varargin{k}; % Set provided axes as place where to plot
    k = k+1; % Set iterator to next input channel
else
    ax = gca; % Set current axes as place where to plot
end

% BE AWARE: *ugly code* below!
% What it does is the following:
% The user can pass agentData to animate() like:
%   animate(data1,data2,date3, ...GLOBAL_PARAMETERS)
% First, the code checks how many agentData have been provided.
% Then, it processes these agentData (more information below).
% Finally, the global parameters are processed.

% extract data from passed arguments
if isa(varargin{k}, 'Simulink.SimulationOutput') % extract simulation data
    simout = varargin{k};
    
    % cell structure to be filled and appended to varargin
    var = {};
    var_opt = {};
    
    % Simulation time
    var{end+1} = simout.tout;

    % target-pose
    if containsData(simout,'gwo')
        var{end+1} = {simout.gwo.signals.values,agentdesign('target','show_trajectory','on')};
    end
    
    % agent-pose & relative target pose
    if containsData(simout,'gwc')
        gwc = simout.gwc.signals.values;
        % estimated target relative pose
        if containsData(simout,'gcobar')
            gcobar = simout.gcobar.signals.values;
            var{end+1} = {gwc,agentdesign('agent','show_trajectory','on'),gcobar,agentdesign('estimate')};
        else
            var{end+1} = {gwc,agentdesign('agent','show_trajectory','on')};
        end
    end
    
    % GP: X data
    if containsData(simout,'X')
        var_opt{end+1} = 'gpdata';
        var_opt{end+1} = {simout.X.time, simout.X.signals.values(:,1:3,:)};
    end
    
    % append to varargin
    varargin = [varargin(1:k) var(:)' varargin(k+1:end) var_opt(:)'];
    
    % set iterator to next input channel
    k = k+1;
end

% check if source time has been given
time_source = [];
if isnumeric(varargin{k}) && isvector(varargin{k})
    time_source = varargin{k};
    k = k+1;
end

% Check how many agents and estimates have been provided
n_rigid = 0; %# agents (+target)
n_rigid_est = 0; %# estimates
k_ = k;
while k_ <= length(varargin) && isAgentData(varargin{k_})
    n_rigid = n_rigid + 1;
    if hasEstimate(varargin{k_})
        n_rigid_est = n_rigid_est + 1;
    end
    k_ = k_+1;
end

% extract data from passed arguments
% BE AWARE: *ugly code* below!
% Every agentData to be processed can have one of the following forms:
%   {gw}, {gw,agentdesign}, {gw,estimate}, {gw,agentdesign,estimate},
%   {gw,estimate,estimatedesign}, {gw,agentdesign,estimate,estimatedesign},
%   {estimate,estimatedesign}, {estimate}
% The order is important!!
% The code below just checks what form the agentData has, and fills two
% cell-arrays for the rigid body motion and the design of every
% agent/estimate.

% body motion of agent/estimate
% form: rigid_gw = {agent1, agent2, ..., agentN, estimate1, estimate 3, ...}
rigid_gw = cell(1,n_rigid+n_rigid_est);
% design of agent/estimate
% form: rigid_design = {agent1, agent2, ..., agentN, estimate1, estimate 3, ...}
rigid_design = cell(1,n_rigid+n_rigid_est);

j_rigid_est = 1; % counter: which estimate is being processed
for i_rigid = 1:n_rigid % loop over all agentData
    agentData = [varargin{k+i_rigid-1} {}]; % enforce a cell layout

    % parse data
    reached_estimate = false;
    l = 1;
    while l <= length(agentData) % if agentdata is of form {gw, ...more}
        if ~reached_estimate % pass pose and design data
            % process agentdata of forms {gw, design, ...} and {gw, ...}
            rigid_gw{i_rigid} = agentData{l};
            if l+1 <= length(agentData) && isa(agentData{l+1},'agentdesign')
                l = l+1;
                rigid_design{i_rigid} = agentData{l};
            else
                rigid_design{i_rigid} = agentdesign();
            end
            reached_estimate = true;
        else % pass pose and design data, but estimate has special case
            % now process all given estimates this agent has
            %   {..., estimate, design, ...}
            rigid_gw{n_rigid+j_rigid_est} = pagemtimes(rigid_gw{i_rigid},agentData{l});
            if l+1 <= length(agentData) && isa(agentData{l+1},'agentdesign')
                l = l+1;
                rigid_design{n_rigid+j_rigid_est} = agentData{l}; % append estimate to end of array
                j_rigid_est = j_rigid_est+1;
            else
                rigid_design{n_rigid+j_rigid_est} = agentdesign('estimate');
            end
        end
        l = l+1;
    end
end

% Process global parameters: {..., parameters}

k = k + n_rigid;

% parse
p = inputParser();
addParameter(p,'gpdata',[],@(x) isGPdata(x))
addParameter(p,'bound',double.empty(3,0),@(x) isnumeric(x) && all(size(x) == [3,2]));
addParameter(p,'projection','perspective',@(x) validatestring(x,{'perspective','orthographic'}));
addParameter(p,'grid','on',@(x) validatestring(x,{'off','on'}));
addParameter(p,'viewangle',3,@(x) isnumeric(x));
addParameter(p,'communication',[],@(x) isnumeric(x) && (isempty(x) || size(x,2) == 2));
addParameter(p,'communication_color','k');
addParameter(p,'communication_style','--',@(x) isa(x,'string'));
addParameter(p,'communication_width',2,@(x) isnumeric(x));
addParameter(p,'fps',20,@(x) isscalar(x));
addParameter(p,'recorder',struct(),@(x) isstruct(x));
parse(p,varargin{k:end});

% process gpdata
time_gp_source = [];
X = [];
if iscell(p.Results.gpdata) && length(p.Results.gpdata) == 2
    time_gp_source = p.Results.gpdata{1};
    X = p.Results.gpdata{2};
end

% process bounds
xbnd = p.Results.bound(1,:);
ybnd = p.Results.bound(2,:);
zbnd = p.Results.bound(3,:);

% process other parameters
projection = p.Results.projection;
showgrid = p.Results.grid;
viewangle = p.Results.viewangle;
communication = p.Results.communication;
communication_color = p.Results.communication_color;
communication_style = p.Results.communication_style;
communication_width = p.Results.communication_width;
fps = p.Results.fps;
recorder_settings = p.Results.recorder;


%% Interpolate time
% if time has been given, ensure smooth playback (even with non-equidistant
% times) of the animation. Solution was found here:
% https://mathworks.com/matlabcentral/answers/414265-interpolate-matrices-for-different-times-in-matlab

% interpolate rigid body motion
time_target = time_source;
if ~isempty(time_source)
    time_target = time_source(1):1/fps:time_source(end);
    for i = 1:n_rigid+n_rigid_est
        rigid_gw{i} = permute(interp1(time_source,permute(rigid_gw{i},[3,2,1]),time_target),[3,2,1]);
    end
end

% interpolate gpdata by holding previous sample point method
if ~isempty(time_gp_source) && ~isempty(X)
    X = permute(interp1(time_gp_source,permute(X,[3,2,1]),time_target,'previous'),[3,2,1]);
end


%% Compute propagation of vertices and patches
% Explanation:
%   The real computation we have to do for a vertical vector w is
%       [w_; 1] =  g * [w; 1]
%   However, the vertices are represented as a matrix of horizontal
%   vectors, thus we calculate instead
%       [w_', 1] = [w', 1] * g'

% vertice datastructures (to be filled)
rigid_vert = cell(1,n_rigid+n_rigid_est);

% propagate rigid body vertices
for i_rigid = 1:n_rigid+n_rigid_est
    rigid_vert{i_rigid} = propagateVertices(rigid_design{i_rigid}.vertices .* rigid_design{i_rigid}.blocksize, rigid_gw{i_rigid});
end

% calculate fitting bounds if not set by user
if isempty(xbnd) && isempty(ybnd) && isempty(zbnd)
    vert_all = cat(1,rigid_vert{:});
    
    % calculate individual bounds
    xbnd = [min(vert_all(:,1,:),[],'all') max(vert_all(:,1,:),[],'all')];
    ybnd = [min(vert_all(:,2,:),[],'all') max(vert_all(:,2,:),[],'all')];
    zbnd = [min(vert_all(:,3,:),[],'all') max(vert_all(:,3,:),[],'all')];
    
    % set bounds to tight scheme
    xbnd = [floor(xbnd(1)) ceil(xbnd(2))];
    ybnd = [floor(ybnd(1)) ceil(ybnd(2))];
    zbnd = [floor(zbnd(1)) ceil(zbnd(2))];
    
    % calculate center
%     all_bound = [min([xbnd(1) ybnd(1) zbnd(1)]), max([xbnd(2) ybnd(2) zbnd(2)])];
%     distance = all_bound(2) - all_bound(1);
%     xcenter = xbnd(2)/2 + xbnd(1)/2;
%     ycenter = ybnd(2)/2 + ybnd(1)/2;
%     zcenter = zbnd(2)/2 + zbnd(1)/2;
    
    % set bounds to make it square
%     xbnd = [floor(xcenter-distance/2), ceil(xcenter+distance/2)];
%     ybnd = [floor(ycenter-distance/2), ceil(ycenter+distance/2)];
%     zbnd = [floor(zcenter-distance/2), ceil(zcenter+distance/2)];
else
    vert_all = cat(1,rigid_vert{:});
    if isempty(xbnd)
        xbnd = [min(vert_all(:,1,:),[],'all') max(vert_all(:,1,:),[],'all')];
    end
    if isempty(ybnd)
        ybnd = [min(vert_all(:,2,:),[],'all') max(vert_all(:,2,:),[],'all')];
    end
    if isempty(zbnd)
        zbnd = [min(vert_all(:,3,:),[],'all') max(vert_all(:,3,:),[],'all')];
    end
end


%% Plotting animation

% plot initial patch
rigid = cell(1,n_rigid);
for i_rigid = 1:n_rigid+n_rigid_est
    linestyle = '-';
    if rigid_design{i_rigid}.edgewidth == 0
        rigid_design{i_rigid}.edgewidth = 1;
        linestyle = 'none';
    end
    rigid{i_rigid} = patch(ax,...
        'Faces',rigid_design{i_rigid}.faces,...
        'Vertices',rigid_vert{i_rigid}(:,:,1),...
        'FaceVertexCData',rigid_design{i_rigid}.color,...
        'FaceColor','interp',...
        'FaceAlpha',rigid_design{i_rigid}.alpha,...
        'LineWidth',rigid_design{i_rigid}.edgewidth,...
        'LineStyle',linestyle);
end

% keep all plots
hold(ax,'on')

% plot gpdata
if ~isempty(X)
    gpdata = plot3(X(:,1,1),X(:,2,1),X(:,3,1),'xr','MarkerSize',10,'LineWidth',2);
end

% plot trajectories
rigid_trajectory = cell(1,n_rigid+n_rigid_est);
for i = 1:n_rigid+n_rigid_est
    if strcmp(rigid_design{i}.show_trajectory,'on')
        trajectory_color = rigid_design{i}.trajectory_color;
        if ~isvector(trajectory_color) % process gradients differently
            rigid_trajectory{i} = patch(...
                'XData',[rigid_gw{i}(1,4,1) nan],...
                'YData',[rigid_gw{i}(2,4,1) nan],...
                'ZData',[rigid_gw{i}(3,4,1) nan],...
                'FaceVertexCData',repeatv(trajectory_color,2),...
                'EdgeColor','interp',...
                'LineWidth',rigid_design{i}.trajectory_width,...
                'LineStyle',rigid_design{i}.trajectory_style);
        else % static color plot
            rigid_trajectory{i} = plot3(rigid_gw{i}(1,4,1),rigid_gw{i}(2,4,1),...
                rigid_gw{i}(3,4,1),...
                'Color',trajectory_color,...
                'LineWidth',rigid_design{i}.trajectory_width,...
                'LineStyle',rigid_design{i}.trajectory_style);
        end
    end
end

% plot communication lines between agents
m_communication = size(communication,1);
rigid_communication = cell(1,m_communication);
for i = 1:m_communication
    agentA = rigid_gw{communication(i,1)};
    agentB = rigid_gw{communication(i,2)};
    rigid_communication{i} = plot3([agentA(1,4,1) agentB(1,4,1)],...
        [agentA(2,4,1) agentB(2,4,1)],[agentA(3,4,1) agentB(3,4,1)],...
        'Color',communication_color,'LineWidth',communication_width,...
        'LineStyle',communication_style);
end

% Show grids on plot
grid(ax,showgrid)

% Limits
xlim(ax,xbnd);
ylim(ax,ybnd);
zlim(ax,zbnd);

% set aspect ratio
daspect(ax,[1 1 1]);

% Set axis-labels
ax.FontSize = 18;
xlabel(ax, 'x [m]', 'FontSize', 20, 'FontWeight', 'bold')
ylabel(ax, 'y [m]', 'FontSize', 20, 'FontWeight', 'bold')
zlabel(ax, 'z [m]', 'FontSize', 20, 'FontWeight', 'bold')

% Show timer
if ~isempty(time_target)
    pos_ = ax.Position;
    t_anno = annotation(ax.Parent,'textbox',[pos_(1)+pos_(3)*0.8 pos_(2)+pos_(4)*0.1 .1 .1],...
        'String',sprintf('t=%.1fs',time_target(1)),'EdgeColor','none',...
        'BackgroundColor',[.5 .5 .5],'FaceAlpha',0.2,'FontSize',13);
end

% Set better 3d projection
set(ax,'Projection',projection)

% View angle
view(ax,viewangle);

% Set 3D-Rotation with mouse ON
rotate3d(ax,'on')

% Draw figure
drawnow

% Set axis as vis3d
% Note: If you do it too early before the figure was rendered, the plot is
%       too close in the front
axis(ax,'vis3d')

% frames to draw
n_frames = size(rigid_vert{1},3);

% video recorder
recorder = [];
if ~isempty(fieldnames(recorder_settings))
    % create path to file if it does not exist yet
    pathparts = strsplit(recorder_settings.file,{'/','\'});
    if length(pathparts) > 1
        filepath = fullfile(pathparts{1:end-1});
        if ~exist(filepath, 'dir'); mkdir(filepath); end
    end
    
    % video recorder settings
    recorder = VideoWriter(recorder_settings.file, recorder_settings.profile);
    recorder.Quality = recorder_settings.quality;
    recorder.FrameRate = fps;
    % ToDo: Insert other video options here
    
    % open video recorder
    fprintf('[animate] recording video to: %s/%s ...', recorder.Path, recorder.Filename);
    open(recorder);
end

% Animation Loop
for k=1:1:n_frames
    t_start = tic;
    
    % update rigid block
    for i_rigid = 1:n_rigid+n_rigid_est
        set(rigid{i_rigid},'Vertices',rigid_vert{i_rigid}(:,:,k));
    end
    
    % update gp data
    if ~isempty(X)
        set(gpdata,'XData',X(:,1,k));
        set(gpdata,'YData',X(:,2,k));
        set(gpdata,'ZData',X(:,3,k));
    end
    
    % update communication lines between agents
    for i = 1:m_communication
        agentA = rigid_gw{communication(i,1)};
        agentB = rigid_gw{communication(i,2)};
        set(rigid_communication{i},'XData',[agentA(1,4,k) agentB(1,4,k)]);
        set(rigid_communication{i},'YData',[agentA(2,4,k) agentB(2,4,k)]);
        set(rigid_communication{i},'ZData',[agentA(3,4,k) agentB(3,4,k)]);
    end
    
    % update trajectory
    for i = 1:n_rigid+n_rigid_est
        if ~isempty(rigid_trajectory{i})
            if isa(rigid_trajectory{i},'matlab.graphics.primitive.Patch') % gradient color trajectories
%                 interval = 1:k;
%                 interval = 1:300:k; % TODO: make this dependent on the distance traveled
                interval = 1:50:k;
                set(rigid_trajectory{i},'XData',[squeeze(rigid_gw{i}(1,4,interval))' nan]);
                set(rigid_trajectory{i},'YData',[squeeze(rigid_gw{i}(2,4,interval))' nan]);
                set(rigid_trajectory{i},'ZData',[squeeze(rigid_gw{i}(3,4,interval))' nan]);
                set(rigid_trajectory{i},'FaceVertexCData',repeatv(rigid_design{i}.color,length(interval)+1));
            else % static color trajectories
                set(rigid_trajectory{i},'XData',rigid_gw{i}(1,4,1:k));
                set(rigid_trajectory{i},'YData',rigid_gw{i}(2,4,1:k));
                set(rigid_trajectory{i},'ZData',rigid_gw{i}(3,4,1:k));
            end
        end
    end
    
    % update timer
    if ~isempty(time_target)
        t_anno.String = sprintf('t=%.1fs',time_target(k));
    end
    
    % Update animation
    if fps > 20
        % 'drawnow' shows every frame as it was written to the figure, but
        % interacting with the figure (rotating, moving, ...) pauses the
        % code. That means: use it for smoother animations e.g. when you
        % want to record a video and every frame is crucial
        drawnow
    else
        % 'drawnow limitrate' only updates the figure in maximal 20fps and
        % only if the renderer is not busy. This option therefore skips a
        % few frames and is not recommended when you record videos, but it
        % has the advantage that you can interact with the figure without
        % the code pausing.
        drawnow limitrate
    end
    
    % append frame to videorecorder
    if ~isempty(fieldnames(recorder_settings))
        % TODO: improve performance with MATLAB backgroundpool as mentioned
        % here: https://mathworks.com/help/matlab/ref/videowriter.html
        writeVideo(recorder, getframe(ax.Parent));
    end
    
    % pause to slow down to desired fps
    t_elapsed = toc(t_start);
    pause(1/fps-t_elapsed);
end

% close video recorder
if ~isempty(recorder)
    close(recorder);
    fprintf(' [done]\n');
end


% -------------- END CODE ---------------
end


function ret = isAgentData(data)
% checks if data is agent data
if iscell(data)
    data = data{1};
end
ret = isnumeric(data) && ndims(data) == 3;
end

function ret = isGPdata(data)
% checks if given data is GP X data (with possibly time attached)

% TODO: Implement me
% isnumeric(x) && any(ndims([]) == [2 3]) && any(size(x,2) == [0 3])
ret = true;
end

function ret = containsData(simdata, field)
fields = simdata.who;
ret = any(matches(fields,field));
end

function vert = propagateVertices(vertices, gw)
% Explanation:
%   The real computation we have to do for a vertical vector w is
%       [w_; 1] =  g * [w; 1]
%   However, the vertices are represented as a matrix of horizontal
%   vectors, thus we calculate instead
%       [w_', 1] = [w', 1] * g'
n_vertices = size(vertices,1); % #vertices
vert = pagemtimes([vertices, ones(n_vertices,1)],'none',gw,'transpose');
vert = vert(:,1:3,:); % Delete 1-column from [vertice, 1]
end

function ret = hasEstimate(agentData)
%%HASESTIMATE checks if agentData provides target estimation data
ret = false;
agentData = [agentData {}]; % enforce a cell layout
for i = 2:length(agentData)
    if isnumeric(agentData{i}) && ndims(agentData{i})==3
        ret = true;
        return;
    end
end
end

function mat = repeatv(mat_,len)
%%REPEATV repeats a matrix vertically starting from the first row again
mat = mat_(mod(0:len-1,size(mat_,1))+1,:);
end