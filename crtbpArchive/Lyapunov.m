% Lyapunov.m
%
% Numerically integrates two test particles in the circular restricted
% three-body problem using the Runge-Kutta-Nystrom 12/10 integrator to
% calculate the Lyapunov exponent.
%
% Program dependences (required to run this code):
%   crtbpRKN1210.m - Integrates 3-body problem
%   rkn1210.m      - This code may be found on the MATLAB file and was
%                    written by Rody P.S. Oldenhuis.
%
% MATLAB-Monkey.com   10/6/2013

clc
clear

nPeriods = 100; % number of orbital periods to run simulation

M1 = 1;        % mass 1
M2 = 0.001;    % mass 2
M = M1 + M2;   % total mass

P = 2*pi * sqrt(1 / M)   % period from Kepler's 3rd law
omega = 2*pi/P;           % angular velocity of massive bodies


times = [0 nPeriods*P];   % set integration limits

R = 1;                    % separation between masses must be 1
r1 = -R*M2/M;             % x coordinate of M1
r2 = R*M1/M;              % x coordinate of M2


LP = LagrangePoints(M2/M);  % get positions of Lagrange points


%%%%%%%%%%  Set initial conditions to match Ceres *** No Chaos
P2 = P*4.60/11.86;      % ratio of Mars's period to Jupiter's
e = 0.079;              % eccentricity
a = R * (P2/P)^(2/3)    % calculate semimajor axis from period
x0 = a*(1-e);           % initial position at perihelion
y0 = 0;                 
vx0 = 0;                % initial velocity at perihelion
vy0 = sqrt(M1*(1+e)/x0);


%%%%%%%%%%  Set initial conditions for 2:3 resonance *** evidence of chaos
% P2 = P*2/3;             % 2:3 resonance
% e = 0.56;               % eccentricity
% a = R * (P2/P)^(2/3);   % calculate semimajor axis from period
% x0 = a*(1+e);           % initial position
% y0 = 0;                 
% vx0 = 0;                % initial velocity
% vy0 = sqrt(M1*(1-e)/x0);

%%%%%%%%%%  Set initial conditions for 7:4 resonance
% P2 = P*7/4;             % 7:4 resonance
% e = 0.1;               % eccentricity
% a = R * (P2/P)^(2/3);   % calculate semimajor axis from period
% x0 = a*(1-e);           % initial position
% y0 = 0;                 
% vx0 = 0;                % initial velocity
% vy0 = sqrt(M1*(1+e)/x0);

%%%%%%%%%%  Set initial conditions 0.025 to right of L4 point
% lp = 4;
% x0 = LP(lp,1) + 0.024;
% y0 = LP(lp,2); 
% v = sqrt(x0^2+y0^2)*omega;
% th = atan2(y0,x0);
% vx0 = v * cos(th+pi/2);
% vy0 = v * sin(th+pi/2);





NP = 2;                 % number of test particles = 2


%%%%%%%%%%  Set plotting flags for integrator
rotatingFlag = true;
animateFlag = false;
trailFlag = true;
PoincareFlag = false;
flags = [rotatingFlag, animateFlag, trailFlag, PoincareFlag];  % plotting flags


%%%%%%%%%%  Integrate 

vals = crtbpRKN1210([M1 M2], [x0 y0 x0 y0-1e-8], [vx0 vy0 vx0 vy0], times, flags);

t = vals(:,1);
pos = vals(:,2:2*NP+1);
vel = vals(:,2*NP+2:4*NP+1);



%%%%%%%%%%  Plot Lyapunov exponents as a function of integration time
figure

subplot(1,2,1)
d0 = sqrt((pos(1,1)-pos(1,3)).^2 + (pos(1,2)-pos(1,4)).^2);
d = sqrt((pos(:,1)-pos(:,3)).^2 + (pos(:,2)-pos(:,4)).^2);
gamma = log(d/d0)./t;
loglog(t,gamma,'b-')
ylabel('\lambda')
xlabel('time')
title(sprintf('Lyapunov Exponent   (m_2/m_1 = %.3f  P/P_0 = %.3f)',M2/M1, P2/P));

subplot(1,2,2)
semilogy(t,d)
xlabel('t')
ylabel('\Delta r')




