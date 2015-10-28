% crtbpRK45.m
%
% Numerically integrate a test particle in the circularly restricted
% three-body problem using the Runge-Kutta 4/5 integrator.
%
% Massive bodies M1 and M2 follow circular orbits around their center of 
% mass point, set to be the origin. The third test particle does affect the
% motion of the two massive bodies. 
%
% To run this code, rkn1210.m must be downloaded from the MATLAB file
% Exchange.  Thanks to Rody P.S. Oldenhuis for his File Exchange code.
%
% We assume G = 1 and R = 1 (distance between the primary bodies)
%
% This function is designed to be called by another function which sets 
% the initial conditions. 
%
% This function calls the following functions:
%   - rkn45.m 
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
%   flag   - 1x1 array containing plotting options (flags are optional):
%            flag(1) = true to create Poincare section
%
% Returned arrays:
%   t      - Tx1 array returning times used in integration, where T =
%            number of time steps
%   pos    - Tx2N array returning (x,y) coordinate pairs for each test
%            particle
%   vel    - Tx2N array returning (vx,vy) velocity component pairs for each 
%            test particle
%   YE     - matrix containing positions and velocities for Poincare
%            Section
%
% MATLAB-Monkey.com   10/18/2013


function [t, pos, vel, YE] = crtbpRK45(masses, pos0, vel0i, times, flag)

if nargin<4 
    display('crtbpRKN1210 requires at least 4 parameters');
    exit;
end

    
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%           Initialization and Parameters               %%%%%%%%%% 
%%%%%%%%%%                                                       %%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%  Flags Controlling Output  %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if nargin >= 5 
    Poincare = flag(1);      % if true, create Poincare section
else
    Poincare = false;
end                 




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


% calculate velocity components in rotating frame
vm = sqrt(vel0i(1)^2 + vel0i(2)^2);
vxu = vel0i(1)/vm;
vyu = vel0i(2)/vm;
vr = omega0 * sqrt(pos0(1)^2+pos0(2)^2);
vnew = vm-vr;
vel0(1) = vxu * vnew;
vel0(2) = vyu * vnew;





%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                      Integrate                        %%%%%%%%%% 
%%%%%%%%%%                                                       %%%%%%%%%%


% Use Runge-Kutta 45 integrator to solve the ODE


if Poincare
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events',@events,'outputfcn',@OutputFcn1);
else
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'outputfcn',@OutputFcn1);
end

% initialize waitbar
wait = waitbar(0, 'integrating...');


% Use Runge-Kutta-Nystrom integrator to solve the ODE
tic;
[t,q,TE,YE] = ode45(@derivatives, times, [pos0 vel0]', options);
DT = toc;

pos = q(:,1:2);
vel = q(:,3:4);




%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                      Functions                        %%%%%%%%%% 
%%%%%%%%%%                                                       %%%%%%%%%%


    %%%%%%%%%%  Derivatives  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Function defining derivatives dx/dt and dy/dt
    % uses the global parameters omega, R10, R20, mu1, mu2 but changeth them not
    function dvdt = derivatives(tf,wf)
           
        R1 = R10;
        R2 = R20;

            dvdt(1) = wf(3);
            dvdt(2) = wf(4);
        
            X = wf(1);            % x component of ii^th particle
            Y = wf(2);              % y component of ii^th particle

            r1 = ((X-R1(1))^2 + (Y-R1(2))^2)^(3/2);
            r2 = ((X-R2(1))^2 + (Y-R2(2))^2)^(3/2);
        
            dvdt(3) = 2*omega0*wf(4) + omega0^2*wf(1) - mu1 * (X-R1(1)) / r1 - mu2 * (X-R2(1)) / r2;
            dvdt(4) = -2*omega0*wf(3) + omega0^2*wf(2) - mu1 * (Y-R1(2)) / r1 - mu2 * (Y-R2(2)) / r2;
        dvdt = dvdt';
    end



    %%%%%%%%%%  Rotation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % rotate column matrix containing (x,y)
    % if dir == 1, rotate forward
    % if dir == -1, rotate backward
    function XYout = rotation(XYf, tf, dir)
        
        s = size(XYf);
        if s(2) > 1
            XYin = XYf';
        else
            XYin = XYf;
        end
        
        tin = tf(end);
        NP = length(XYin)/2;
        rot = [cos(omega0*tin) -sin(omega0*tin);
               sin(omega0*tin)  cos(omega0*tin)];

        if dir == -1
            rot = rot';
        end
        
        for ii = 1:NP
            XYout = (rot * XYin)';
        end
    end


    %%%%%%%%%%  Events - Poincare Section  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %

    function [value,isterminal,direction] = events(tin,y,dy)
        % Locate the time when y1 passes through zero in an
        % increasing direction

        value = y(2);     % Detect when particle passes y1 = 0
        isterminal = 0;   % Don't stop the integration
        direction = -1;        % Positie direction only
    end



    %%%%%%%%%%  OutputFcn1:  Status Bar  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % the output function
    function stop = OutputFcn1(t,y,flag)%#ok
        % don't stop
        stop = false;
        % only after sucessfull steps
        if ~isempty(flag), return
        else
            wait = waitbar(t(end)/times(end), wait);
        end
    end


end
