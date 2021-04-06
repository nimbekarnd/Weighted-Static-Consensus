%% Weighted Static Consensus over Line, Cycle, and Complete Graphs
% This script executes the consensus protocol over N agents with weights
% intended to ensure the agents maintain a separation distance from each
% other. The network topology can be set to that of a line, cycle or
% complete graph.
%
% Author: Yancy Diaz-Mercado
% Date Created: 2020/10/6
%
% Copyright (c) 2020
% Collaborative Controls and Robotics Laboratory
% University of Maryland, College Park
%
% All rights reserved.
%
%% Initialize Robotarium Object
% Before running this script, you must have downloaded the Robotarium API
% and run the init.m script contained with it. This will add the Robotarium
% class and functions to the MATLAB path. This must be done whenever the
% path is reset (e.g., when MATLAB is restarted).
% 
% Number of Agents
numAgents = 8;
% Node desired separation
desiredSep = 0.5;
% Select a graph topology. 
graphTopology = 'line'; % Choose between: 'line', 'cycle', 'complete'
% Initial configuration for the agents
randOrient = linspace(0,2*pi*(numAgents-1)/numAgents,numAgents);
randDistX = 1.3 + 0.2*rand(1,numAgents);
randDistY = 0.7 + 0.2*rand(1,numAgents);
initialPose = [randDistX.*cos(randOrient);randDistY.*sin(randOrient);randOrient-pi];
% Attempt to instatiate a Robotarium object
try
    r = Robotarium('NumberOfRobots', numAgents, 'ShowFigure',true,...
        'InitialConditions',initialPose);
catch e
    if strcmp(e.identifier,'MATLAB:UndefinedFunction')
        warning('MATLAB was unable to find the Robotarium API. Did you run the init.m script prior to execution?')
    end
    e.throw
end
%% Visualization
% Get figure handle
figureHandle = r.figure_handle;
hold on
% Relabel node connectivity based on orientation to reduce crossing edges
x0 = r.get_poses;
switch graphTopology
    case 'line'
        % Define the line graph topology
        A = diag(ones(numAgents-1,1),1) + diag(ones(numAgents-1,1),-1);
    case 'cycle'
        % Define the cycle graph topology
        A = diag(ones(numAgents-1,1),1) + diag(ones(numAgents-1,1),-1);
        A(1,numAgents) = 1;
        A(numAgents,1) = 1;
    case 'complete'
        % Define the complete graph topology
        A = ones(numAgents) - eye(numAgents);
    otherwise
        error('Must select a graph topology in these three classes: ''line'', ''cycle'', or ''complete''.')
end
agentInd = 1:numAgents;
% Get edge indices for plotting dynamically
numEdges = sum(A,'all')/2;
edgeIndices = zeros(numEdges,2);
edgeCount = 1;
for ii = 1:numAgents-1
    for jj = ii+1:numAgents
        if A(ii,jj) == 1
            edgeIndices(edgeCount,1) = ii;
            edgeIndices(edgeCount,2) = jj;
            edgeCount = edgeCount + 1;
        end
    end
end
% Handle for agent interaction edges
edgeX = [x0(1,edgeIndices(:,1));x0(1,edgeIndices(:,2));nan(1,numEdges)];
edgeY = [x0(2,edgeIndices(:,1));x0(2,edgeIndices(:,2));nan(1,numEdges)];
edgeHandle = plot(edgeX(:),edgeY(:),'k','LineWidth',3);
% Linear parameterization for plotting
linPar = linspace(0,1);
% Agent size
agentSize = r.robot_diameter;
% Pre-allocate handles
agentHandle = gobjects(numAgents,1);
agentLabelHandle = gobjects(numAgents,1);
for ii = 1:numAgents
    agentHandle(ii) = patch(...
        'XData',x0(1,ii)+agentSize*cos(2*pi*linPar),...
        'YData',x0(2,ii)+agentSize*sin(2*pi*linPar),...
        'FaceColor',0.25*[1 1 1],'EdgeColor','k','LineWidth',2,'FaceAlpha',0.25);
    agentLabelHandle(ii) = text(...
        x0(1,ii),x0(2,ii),num2str(ii),'Color','w',...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
end
r.step;
%% Control Utilities
% Single integrator to unicycle transformation
si2uni = create_si_to_uni_mapping;
% Barrier certificate
barrierControl = create_si_barrier_certificate('SafetyRadius',1.25*agentSize);
% Max speed saturation
maxSpeedLim = 0.25*r.max_linear_velocity;
%% Run Weighted Consensus
dx = ones(2,numAgents);
% Main loop
while any(vecnorm(dx) > 5e-3)
    % Update pose information
    x = r.get_poses;
    for ii = 1:numAgents
        % Update agent visualization
        set(agentHandle(ii),...
            'XData',x(1,ii)+agentSize*cos(2*pi*linPar),...
            'YData',x(2,ii)+agentSize*sin(2*pi*linPar));
        set(agentLabelHandle(ii),'Position',[x(1,ii),x(2,ii)]);
    end
    % Update edge visualization
    edgeX = [x(1,edgeIndices(:,1));x(1,edgeIndices(:,2));nan(1,numEdges)];
    edgeY = [x(2,edgeIndices(:,1));x(2,edgeIndices(:,2));nan(1,numEdges)];
    set(edgeHandle,'XData',edgeX(:),'YData',edgeY(:))
    % Render
    drawnow limitrate
    % Reset velocity vector
    dx = zeros(2,numAgents);
    %%% Compute the control
    for ii = 1:numAgents
        % Own state
        xi = x(1:2,ii);
        for jj = agentInd(A(ii,:)==1)
            % Neighbor state
            xj = x(1:2,jj);
            % Distance to neighbor
            distij = norm(xi-xj);
            % Inter-agent weight
            wij = distij-desiredSep; % <-- Change this value to achieve an inter-agent spacing 'distij' of 'desiredSep'
            % Compute control law
            dx(:,ii) = dx(:,ii) + wij*(xj-xi);
        end
    end    
    % Wrap controller in barrier certificate for safety
    dx =  barrierControl(dx,x);
    % Convert the controller to unicycle
    dxu = si2uni(dx,x);
    % Scale speed
    maxSpeed = max(abs(dxu(1,:)));
    if maxSpeed > maxSpeedLim
        dxu = dxu*maxSpeedLim/maxSpeed;
    end
    % Send velocity commands
    r.set_velocities(agentInd,dxu);
    r.step
end
% We're done!
disp('Agents have converged!')
text(0,-0.75,{'Agents have converged!';['Min Neighbor Dist = ',num2str(min(A.*sqrt((x(1,:)-x(1,:).').^2+(x(2,:)-x(2,:).').^2)+10*(A==0),[],'all'))]},'HorizontalAlignment','center','FontSize',18)
pause(3)
r.debug