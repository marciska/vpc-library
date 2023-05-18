classdef agentdesign
    %AGENTDESIGN options for agents used by animate().
    %   Sets default values for agent design like Verices, faces, colors
    %   that can be altered.
    %
    %   Property of: Fujita-Yamauchi Lab, University of Tokyo, 2021
    %   e-mail: marcoomainska@g.ecc.u-tokyo.ac.jp
    %   Website: https://www.scl.ipc.i.u-tokyo.ac.jp
    %   December 2021
    
    properties
        blocksize (1,3) double {mustBePositive, mustBeFinite} = .4*ones(1,3)
        vertices (:,3) double {mustBeReal, mustBeFinite} = [-.5  -.3  -.5;
                                                             .5  -.3  -.5;
                                                             .5  -.3   .5;
                                                            -.5  -.3   .5;
                                                              0   .7    0]
        alpha (1,1) double {mustBeInRange(alpha,0,1)} = 1;
        edgewidth (1,1) double {mustBeNonnegative, mustBeFinite} = 1;
        show_trajectory {mustBeMember(show_trajectory,{'off','on'})} = 'off';
        trajectory_color = "inherit";
        trajectory_style string = '-';
        trajectory_width double {mustBePositive} = 2;
    end
    properties(Dependent)
        faces
        color
    end
    properties(Access=private)
        faces_ (:,:) double = [1 2 3 4;
                              1 2 5 nan;
                              2 3 5 nan;
                              3 4 5 nan;
                              4 1 5 nan]
        color_ (:,3) double {mustBeNonnegative} = ones(5,1).*[0 0.4470 0.7410];
    end
    
    methods
        function this = agentdesign(varargin)
            %AGENTDESIGN constructs an instance of this class
            %   Detailed explanation goes here
            %obj.Property1 = inputArg1 + inputArg2;
            
            % Process scheme on first argument
            scheme = 'default';
            start_idx = 1;
            expected_schemes = {'default', 'target', 'estimate', 'agent', 'sphere'};
            if nargin >=1 && any(strcmp(varargin{1},expected_schemes))
                scheme = varargin{1};
                start_idx = start_idx + 1;
            end
            switch scheme
                case 'target'
                    this.color = hsv(size(this.vertices,1));
                    this.edgewidth = 0;
                case 'estimate'
                    this.color = ones(size(this.vertices,1),1).*[0.4660 0.6740 0.1880];
                    this.alpha = .5;
                case 'agent'
                    this.color = ones(size(this.vertices,1),1).*[0 0.4470 0.7410];
                case 'sphere'
                    [x, y, z] = sphere(20);
                    f = figure('Visible', 'off');
                    ax = axes(f);
                    h = surf(ax,x,y,z);
                    [faces,vertices,~] = surf2patch(h);
                    close(f);
                    this.vertices = vertices;
                    this.faces = faces;
                    % this.color = hsv(size(this.vertices,1));
                    this.color = ones(size(this.vertices,1),1).*[0.8863 0.5843 0.4706];
                    this.edgewidth = 0;
                otherwise
                    this.color = [0 0.4470 0.7410];
            end
            
            % Parse
            p = inputParser();
            addParameter(p,'vertices',this.vertices);
            addParameter(p,'faces',this.faces);
            addParameter(p,'blocksize',this.blocksize);
            addParameter(p,'color',this.color);
            addParameter(p,'alpha',this.alpha);
            addParameter(p,'edgewidth',this.edgewidth);
            addParameter(p,'show_trajectory',this.show_trajectory);
            addParameter(p,'trajectory_color',this.trajectory_color);
            addParameter(p,'trajectory_style',this.trajectory_style);
            addParameter(p,'trajectory_width',this.trajectory_width);
            parse(p,varargin{start_idx:end})
            
            % Set properties
            this.vertices = p.Results.vertices;
            this.faces = p.Results.faces;
            this.blocksize = p.Results.blocksize;
            switch scheme
                case 'target'
                    this.color = hsv(size(this.vertices,1));
                case 'estimate'
                    this.color = ones(size(this.vertices,1),1).*[0.4660 0.6740 0.1880];
                    this.alpha = .5;
                case 'agent'
                    this.color = ones(size(this.vertices,1),1).*[0 0.4470 0.7410];
                otherwise
                    this.color = [0 0.4470 0.7410];
            end
            if any(size(p.Results.color,1) == [1 size(this.vertices,1)])
                if size(p.Results.color,1) == 1
                    this.color = ones(size(this.vertices,1),1).*p.Results.color;
                else
                    this.color = p.Results.color;
                end
            end
            this.alpha = p.Results.alpha;
            this.edgewidth = p.Results.edgewidth;
            this.show_trajectory = p.Results.show_trajectory;
            this.trajectory_color = p.Results.trajectory_color;
            this.trajectory_style = p.Results.trajectory_style;
            this.trajectory_width = p.Results.trajectory_width;
        end
        
        function this = set.vertices(this,vertices)
            this.vertices = vertices;
        end
        
        function this = set.faces(this,faces)
            assert(max(faces(:)) <= size(this.vertices,1),'Mismatch in vertices and faces! First specify vertices, then the matching faces.');
            this.faces_ = faces;
        end
        function faces = get.faces(this)
            faces = this.faces_;
        end
        
        function this = set.color(this,color)
            if size(color,1) == 1
                color = repmat(color, size(this.vertices,1), 1);
            end
            assert(size(color,1) == size(this.vertices,1),'Mismatch in vertices and colors! First specify vertices, then the matching color.');
            this.color_ = color;
        end
        function color = get.color(this)
            color = this.color_;
        end
        
        function this = set.blocksize(this,blocksize)
            this.blocksize = blocksize;
        end
        
        function this = set.alpha(this,alpha)
            this.alpha = alpha;
        end

        function this = set.edgewidth(this,edgewidth)
            this.edgewidth = edgewidth;
        end
        
        function this = set.show_trajectory(this,show_trajectory)
            this.show_trajectory = show_trajectory;
        end
        
        function this = set.trajectory_color(this,trajectory_color)
            this.trajectory_color = trajectory_color;
        end
        function trajectory_color = get.trajectory_color(this)
            trajectory_color = this.trajectory_color;
            % check if color must be inherited from object
            if strcmp(trajectory_color,'inherit')
                trajectory_color = mean(this.color);
            elseif strcmp(trajectory_color,'inherit-gradient')
                trajectory_color = this.color;
            end
            % check if color is a matrix of color codes
            if ~isvector(trajectory_color)
                if nnz(diff(trajectory_color,1))==0 % all rows are equal -> same color for entire object
                    trajectory_color = trajectory_color(1,:);
                end
            end
        end
        
        function this = set.trajectory_style(this,trajectory_style)
            this.trajectory_style = trajectory_style;
        end
        
        function this = set.trajectory_width(this,trajectory_width)
            this.trajectory_width = trajectory_width;
        end
    end
end

