% crtbRKN1210.m
%
% Numerically integrate a test particle in the circularly restricted
% three-body problem using the Runge-Kutta-Nystrom 12/10 integrator.
%
% Massive bodies M1 and M2 follow circular orbits around their center of 
% mass point, set to be the origin. The third test particle does affect the
% motion of the two massive bodies. 
%
% We assume G = 1 and R = 1 (distance between the primary bodies)
%
% This function is designed to be called by another function which sets 
% the initial conditions. 
%
% This function calls the following functions:
%   - rkn1210.m (written by Rody P.S. Oldenhuis, MATLAB File Exchange)
%   - LagrangePoints.m
%
% Passed parameters:
%   masses - 1x2 matrix containing m1 and m2
%   pos0   - 1x2N matrix cointaining pairs of (x,y) coordinates for
%            the initial conditions of the test particles. If one test
%            particle is used, pos0 will be a 1x2 matrix containing (x0, y0)
%            For 2 test particles, posq will contain (x10, y10, x20, y20)...
%   vel0   - 1x2N matrix containing pairs of (vx, vy) values for each of
%            the initial velocities of the test particles. 
%   times  - 1x2 array containing begin and end times [tBegin, tEnd]
%   flag   - 1x3 array containing plotting options (flags are optional):
%            flag(1) = true to plot and return values in rotating frame
%            flag(2) = true to animate the orbit as it integrates. Cool to
%                      watch, but it slows the calculate way down.
%            flag(3) = true to leave a trail behind the test particles when
%                      animation is turned on.
%
% Returned arrays:
%   t      - Tx1 array returning times used in integration, where T =
%            number of time steps
%   pos    - Tx2N array returning (x,y) coordinate pairs for each test
%            particle
%   vel    - Tx2N array returning (vx,vy) velocity component pairs for each 
%            test particle
%
%
% MATLAB-Monkey.com   10/6/2013


function returnVals = crtbpRKN1210(masses, pos0, vel0, times, flag)

if nargin<4 
    display('crtbpRKN1210 requires at least 4 parameters');
    exit;
end

    
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%           Initialization and Parameters               %%%%%%%%%% 
%%%%%%%%%%                                                       %%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%  Flags Controlling Output  %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if nargin >= 5
    rotatingFrame = flag(1); % if true, transform to rotating frame
                             % if false, display inertial frame

    animate = flag(2);       % if true, animate motion
                             % if false, just display final trajectory

    leaveTrail = flag(3);    % if true and animation is on, past positions will
                             % not be erased.  If false and animation is on,
                             % no 'trail' of past positions will be left 
end                 

progressBar = false;     % if true, show progress bar

animateDelay = 0.0;      % delay (in seconds) between frames of animation.
                         % smaller the value, faster the animation will play
                        
leaveTrail = flag(3);    % if true and animation is on, past positions will
                         % not be erased.  If false and animation is on,
                         % no 'trail' of past positions will be left 
                        
trailLength = 700;       % if leaveTrail is true, this sets the length of
                         % trail




%%%%%%%%%%%%%%%%%%%%%%%%%  System Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NP = length(pos0)/2;  % number of test particles.  Assumes 2D

M1 = masses(1);  % mass 1
M2 = masses(2);  % mass 2
M = M1 + M2;     % total mass

G = 1;           % Gravitational Constant
R = 1;           % distance between M1 and M2 set to 1



%%%%%%%%%%%%%%%%  Orbital properties of two massive bodies  %%%%%%%%%%%%%%%
mu = G*M;
mu1 = G*M1;
mu2 = G*M2;

R10 = [-M2/M; 0];           % initial position of M1
R20 = [M1/M; 0];            % initial position of M2

P = 2*pi * sqrt(R^3 / mu);  % period from Kepler's 3rd law
omega0 = 2*pi/P;            % angular velocity of massive bodies

plotLimits = 1.5;           % limit for plotting
tEnd = times(end);


LP = LagrangePoints(M2/M);   % calculate Lagrange points




%%%%%%%%%%%%%%%%%%%%%%%%%%%   Trail Properties  %%%%%%%%%%%%%%%%%%%%%%%%%%%

trail = ones(trailLength,2*NP);
for j = 1:2*NP
    trail(:,j) = pos0(j)*ones(trailLength,1);
end


%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                      Integrate                        %%%%%%%%%% 
%%%%%%%%%%                                                       %%%%%%%%%%


% Use Runge-Kutta 45 integrator to solve the ODE

if animate
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'outputfcn',@OutputFcn2);
else
    % initialize waitbar
    if progressBar
        wait = waitbar(0, 'integrating...');
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'outputfcn',@OutputFcn1);
    else
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    end
end

% Use Runge-Kutta-Nystrom integrator to solve the ODE
tic;
% [t,pos,vel] = rkn1210(@derivatives, times, pos0, vel0, options);
[t,pos,vel,TE,YE] = rkn1210(@derivatives, times, pos0, vel0, options);

DT = toc;

if rotatingFrame
    for j=1:length(t)
        pos(j,:) = rotation(pos(j,:)',t(j),-1);
    end
end


% returnVals = [t, pos, vel, YE];
returnVals = [t, pos, vel];




%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                      Functions                        %%%%%%%%%% 
%%%%%%%%%%                                                       %%%%%%%%%%


    %%%%%%%%%%  Derivatives  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Function defining derivatives dx/dt and dy/dt
    % uses the global parameters omega, R10, R20, mu1, mu2 but changeth them not
    function dvdt = derivatives(tf,wf)
        NP = length(wf)/2;      % number of particles (assumes 2D)     
                
        rot = [cos(omega0*tf) -sin(omega0*tf);
               sin(omega0*tf)  cos(omega0*tf)];
           
        R1 = rot * R10;
        R2 = rot * R20;

           
        for ii = 1:NP
            X = wf(2*ii-1);            % x component of ii^th particle
            Y = wf(2*ii);              % y component of ii^th particle

            r1 = ((X-R1(1))^2 + (Y-R1(2))^2)^(3/2);
            r2 = ((X-R2(1))^2 + (Y-R2(2))^2)^(3/2);
        
            dvdt(2*ii-1) = - mu1 * (X-R1(1)) / r1 - mu2 * (X-R2(1)) / r2;
            dvdt(2*ii) =  - mu1 * (Y-R1(2)) / r1 - mu2 * (Y-R2(2)) / r2;
        end
    end



    %%%%%%%%%%  OutputFcn1:  Status Bar  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % the output function
    function stop = OutputFcn1(t,y,dy,flag)%#ok
        % don't stop
        stop = false;
        % only after sucessfull steps
        if ~isempty(flag), return
        else
            wait = waitbar(t/tEnd, wait);
        end
    end



    %%%%%%%%%%  OutputFcn2:  Animation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % the output function for making the animation
    function stop = OutputFcn2(t,y,dy,flag)
        % don't stop
        stop = false;
        % only after sucessfull steps
        if ~isempty(flag), return
        
        else
            NP = length(y)/2;      % number of particles (assumes 2D)     

            if rotatingFrame
                
                XYrot = rotation(y,t,-1);
                
                if leaveTrail
                    for ii = 1:NP
                        trail(:,2*ii-1) = [XYrot(2*ii-1); trail(1:end-1,2*ii-1)];
                        trail(:,2*ii) =   [XYrot(2*ii); trail(1:end-1,2*ii)];

                        plot(trail(:,2*ii-1),trail(:,2*ii),'-');
                        hold on
                    end
                end
                
                for ii = 1:NP
                    plot(XYrot(2*ii-1), XYrot(2*ii), '.')
                    hold on
                end
                plot(LP(:,1),LP(:,2),'k+');
                plot(R10(1),R10(2),'ro','MarkerFaceColor','r');
                plot(R20(1),R20(2),'go','MarkerFaceColor','g');
            else
                if leaveTrail
                    for ii = 1:NP
                        trail(:,2*ii-1) = [y(2*ii-1); trail(1:end-1,2*ii-1)];
                        trail(:,2*ii) = [y(2*ii); trail(1:end-1,2*ii)];

                        plot(trail(:,2*ii-1),trail(:,2*ii),'-');
                        hold on
                    end
                end
                for ii = 1:NP
                    plot(y(2*ii-1), y(2*ii), '.')
                    hold on
                end
                R1 = rotation(R10,t,1);
                R2 = rotation(R20,t,1);
                
                plot(R1(1),R1(2),'ro','MarkerFaceColor','r');
                plot(R2(1),R2(2),'go','MarkerFaceColor','g');
            end
            
            axis equal
            axis(plotLimits*[-1 1 -1 1])
            hold off
            pause(animateDelay), drawnow % flush plotting commands
        end
    end



    %%%%%%%%%%  Rotation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % rotate column matrix containing (x,y)
    % if dir == 1, rotate forward
    % if dir == -1, rotate backward
    function XYout = rotation(XYin, tin, dir)
        
        NP = length(XYin)/2;
        rot = [cos(omega0*tin) -sin(omega0*tin);
               sin(omega0*tin)  cos(omega0*tin)];
        
        if dir == -1
            rot = rot';
        end
        
        for ii = 1:NP
            XYout(2*ii-1:2*ii) = (rot * XYin(2*ii-1:2*ii))';
        end
    end
end
