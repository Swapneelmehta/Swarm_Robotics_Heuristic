function [P,Nit] = Project(N,NIter,Bound,Alpha,Beta,R,proximity)
% Default Vaules

if nargin < 7
    proximity=2;
    if nargin < 6
        R(1)=100;
        R(2)=2;
        R(3)=5;
        if nargin < 5
            % Local search coefficient
            Beta(1) = 2.0; ...beta1
                % Social or Global search coefficient
            Beta(2) = 2.0; ...beta2
                % Diversity Factor
            Beta(3) = 1.0;....gamma
             if nargin < 4
                % The inertial coefficient
                Alpha(1) = 1; ...alpha
                % Daming in inertia
                Alpha(2)  = 1; ....alpha_damp
                if nargin < 3
                    % Upper Bound
                    Bound(1) = 100;
                    % Lower Bound
                    Bound(2) = -100;
                    U=Bound(1);
                    L=Bound(2);
                    if nargin < 2
                        % Maximum number of Iterations
                        NIter = 300;
                        if nargin < 1
                        % Population size of the swarm particles
                        N = 10;
                        end
                    end
                end
            end
        end
    end
end


%% User define parameters
U = Bound(1);        % Upper Bound
L = Bound(2);        % Lower Bound
alpha = Alpha(1);       % The inertial coefficient
alpha_damp  = Alpha(2); % Daming in intertia
beta1 = Beta(1);        % Local search coefficient
beta2 = Beta(2);        % Social or Global search coefficient
gamma = Beta(3);        % Diversity Factor
V = zeros(N,2);         % Velocity of the particles
Vmax=1;         % Scalling facotor of the Velocity
r1=R(2);
r2=R(1);
r3=R(3);        %
r4=proximity;   % proximity of the particle to detect the object

%% Set the position of the initial swarm

X= rand(N,2);
P = (U-L)*(rand(1,2)-0.5);
X1(:,1) = 10.*X(:,1) ; % Scaled Population
X1(:,2) = 10.*X(:,2) + L;
Y = zeros(N,1);
Y(:,1) = inf;
PL = (U-L).*lhsdesign(N,2)+L;              % personal best value
tic;
P2=[];
obstacle=[];
History=[];

%% Walls
walls=[L,L,L,U;L,U,U,U;U,U,U,L;U,L,0.2*U,L;
    0.2*L,L,L,L;U,0,0.2*U,0;0.2*L,0,L,0;
    U,0.5*L,0.2*U,0.5*L;0.2*L,0.5*L,L,0.5*L;
    U,0.5*U,0.2*U,0.5*U;0.2*L,0.5*U,L,0.5*U;
    0.2*L,0.75*L,0.2*L,0.5*L;0.2*L,0,0.2*L,0.25*L;
    0.2*L,0.5*U,0.2*L,0.25*U;0.2*L,U,0.2*L,0.75*U;
    0.2*U,0.75*L,0.2*U,0.5*L;0.2*U,0,0.2*U,0.25*L;
    0.2*U,0.5*U,0.2*U,0.25*U;0.2*U,U,0.2*U,0.75*U];

B = [L,L;U,U;L,U;U,L;
    0,0;0,L;0,U;L,0;U,0;
    L,0.5*U;0,0.5*U;U,0.5*U;
    L,0.5*L;0,0.5*L;U,0.5*L];

for i1=1:length(walls)
    P1 = A_wall(walls(i1,:));
    P2 = [P2;P1];
end

dist=bsxfun(@hypot,B(:,1)-P(1),B(:,2)-P(2));
out = B(dist==min(dist),:);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% The main PSO Algorithm%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for Nit=1:NIter
    History=[History X1];
    % Proximity form the wall of each partical
    for i3 =1:N
        dist=bsxfun(@hypot,P2(:,1)-X1(i3,1),P2(:,2)-X1(i3,2));
        proximity(i3,1) = min(dist);
        P21 = P2(dist==min(dist),:);
        CollisionPoint(i3,:)=P21(1,:);
        if proximity(i3,1)<r4
            obstacle=[obstacle;P21];
            obstacle = unique(obstacle,'rows');
        end
    end
    
    % Evaluate the function
    for i3 = 1:N
        fun_val(i3,1) = lin(X1(i3,:),P);
    end
    
    % Compare the function values to find the best ones
    for i4 = 1:N
        if fun_val(i4,1) <= Y(i4)
            if fun_val(i4,1)<r2
                PL(i4,:) = X1(i4,:);      % update personal best  position,
            end
            Y(i4) = fun_val(i4,1);    % update the best value so far
        end
    end
    
    out = find(fun_val > r2);
    Y1 =sortrows(Y);
    
    % find the global best function value
    [PG,pg] = min(Y);
    DisM=Dis(X1);
    Y1=sum(DisM(:,:));
    [~,pg1] = min(Y1);
    
    % update the velocity of the particles
    if length(out) < N   % if any single partticle enter the range of the objective
        disv1= (rand(N,1)*gamma).*(X1(pg,1)-X1(:,1));
        disv2= (rand(N,1)*gamma).*(X1(pg,2)-X1(:,2));
        V(:,1) = alpha*V(:,1)+(rand(N,1)*beta1).*(PL(:,1)-X1(:,1))+(rand(N,1)*beta2).*(X1(pg,1)-X1(:,1))- rand(N,1).*disv1;
        V(:,2) = alpha*V(:,2)+(rand(N,1)*beta1).*(PL(:,2)-X1(:,2))+(rand(N,1)*beta2).*(X1(pg,2)-X1(:,2))- rand(N,1).*disv2;
    else  % If the particles have not in the range of the objective (No Global LeadAer)
        disv1= (rand(N,1)*gamma).*(X1(pg1,1)-X1(:,1));
        disv2= (rand(N,1)*gamma).*(X1(pg1,1)-X1(:,1));
        V(:,1) = alpha*V(:,1)+(rand(N,1)*beta1).*(PL(:,1)-X1(:,1))-disv1;
        V(:,2) = alpha*V(:,2)+(rand(N,1)*beta1).*(PL(:,2)-X1(:,2))-disv2;
    end
    
    % Normalizing The Velocity Vector
    for i=1:N
        w(i,1)=sqrt(V(i,1)^2+V(i,2)^2) ;
        if w(i,1)==0
            w(i,1)=1;
        end
        V(i,:)=Vmax.*V(i,:)/w(i);
    end
    
    % Deteact the collision points if the obstacle is in the proximity of any particle
    for i4 = 1:N
        if proximity(i4,1)< r4
            u = V(i4,:);
            v = X1(i4,:)-CollisionPoint(i4,:);
            angle=acosd(dot(u,v)/(norm(u)*norm(v)));
            if angle <=45
                X1(i4,:) = X1(i4,:)-(rand(1)+1)*V(i4,:)+ 4*(X1(i4,:)-CollisionPoint(i4,:))/sqrt(proximity(i4,1));
            elseif angle <90 && angle > 45
                X1(i4,:) = X1(i4,:)-(rand(1)+1)*V(i4,:)+ 4*(X1(i4,:)-CollisionPoint(i4,:))/sqrt(proximity(i4,1));
            elseif angle <160 && angle > 90
                X1(i4,:) = X1(i4,:)-2*V(i4,:)+ 4*(X1(i4,:)-CollisionPoint(i4,:))/sqrt(proximity(i4,1));
            else
                X1(i4,:) = X1(i4,:)+ 2*rand(1)*(X1(i4,:)-CollisionPoint(i4,:));
            end
        end
    end
    
    
    % Plot the graph
    clf;
    % Plot the wall boundaries in the graph
    for i=1:length(walls)
        wall(walls(i,:))
        hold on;
    end
    % PLot the detected portion of the wall by the particles
    if length(obstacle)>1
        plot(obstacle(:,1),obstacle(:,2),'kx')
    end
    % Plot the particle motion in each Iteration
    plotfun(X1,P,Nit,U,r1,r2,r3)
    axis([L U L U])
    pause(0.05)
    
    % update the position by the velocity each particle
    for i2 =1:2
        X1(:,i2) = X1(:,i2) + V(:,i2);
    end
    % To privent the explosion of the swarm
    % Population is contraint to be in Upper and lower bound
    UP = find(X1 > U);
    for i7 = UP
        X1(i7) = U;
    end
    LP = find(X1 < L);
    for i8 = LP
        X1(i8) = L;
    end
    
    
    
    % Update in the inertia alpha
    if PG < r3
        alpha=alpha*alpha_damp;
    end
    
    % Termination Criteria to stop the algorithm
    if  PG < r1
        break;
    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Quadratic function for the gas leak

function [Y]=lin(xx,P)
d=length(xx);
sum=0;
for ii = 1:d
    xi = xx(ii)-P(ii);
    sum = sum + (xi^2);
end
Y = sqrt(sum);
end

% funciton to calculate the distance between every particle
function DisM=Dis(x)
N=length(x);
for i=1:N
    for j=1:N
        DisM(i,j)=sqrt((x(i,1)-x(j,1))^2+(x(i,2)-x(j,2))^2);
    end
end
end

% Plot the figures
function plotfun(X1,P,Nit,U,r1,r2,r3)
L= -U;
xv = [U;L;L;U;U];
yv = [U;U;L;L;U];
th = 0:pi/100:2*pi;
x1 = r1 * cos(th) + P(1);
y1 = r1 * sin(th) + P(2);
[in1,~] = inpolygon(x1,y1,xv,yv);
x2 = r2 * cos(th) + P(1);
y2 = r2 * sin(th) + P(2);
[in2,~] = inpolygon(x2,y2,xv,yv);
th1 = 0:pi/10:2*pi;
x3 = r3 * cos(th1) + P(1);
y3 = r3 * sin(th1) + P(2);
hold on;
plot(x1(in1), y1(in1),'r');
plot(x2(in2), y2(in2),'.');
plot(x3, y3,'.');
plot(X1(:,1),X1(:,2), 'bx');
plot(P(1),P(2), 'rx');
title(['Itretion No.',num2str(Nit)])
xlabel('Y')
ylabel('X')
end

% Just to draw the wall in the figure
function wall(X0)
X11 = [X0(1), X0(3)];
Y11 = [X0(2), X0(4)];
plot(X11,Y11,'b','Linewidth',1)
end

% Discretize the wall into large number of point which will be
% the obstacle
function [P1]=A_wall(X0)
X1 = linspace(X0(1), X0(3),100)';
Y1 = linspace(X0(2), X0(4),100)';
P1 = [X1 Y1];
end