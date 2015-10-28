% LagrangePlot.m
%
% Plots Lagrange points and zero-velocity curves that pass through them
% for the circular, restricted 3-body problem. 
%
% Assumes G = 1 and r = 1, where r  is the distance between masses m1 & m2.
% The origin of the coordinate system is placed at the center-of-mass point
%
% Program dependencies:
%    crtbpPotential.m - returns pseudo-potential
%    lagrangePoints.m - returns 5x3 array containing (x,y,z) coordinates
%                       of the Lagrange points for given values of m1, m2
%
% Matlab-Monkey.com  10/10/2013



%%%%%%%%%%%%%%%%%%%%%  Parameters and Initialization  %%%%%%%%%%%%%%%%%%%%%

M1 = 1;      % mass 1
M2 = 0.1;    % mass 2
M = M1 + M2; % total mass

P = 2*pi * sqrt(R^3 / mu);  % period from Kepler's 3rd law
omega0 = 2*pi/P;            % angular velocity of massive bodies

% find Lagrange points and the pseudo-potential
LP = lagrangePoints(M2/(M1+M2))
LP1_level = crtbpPotential(mu1, mu2, LP(1,1), LP(1,2));
LP2_level = crtbpPotential(mu1, mu2, LP(2,1), LP(2,2));
LP3_level = crtbpPotential(mu1, mu2, LP(3,1), LP(3,2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%     Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot zero-velocity curves that run through L1, L2, L3
contour(X,Y,U,[LP1_level LP2_level LP3_level]);
hold on

% plot m1 (yellow dot) and m2 (green dot)
plot(-M2/(M1+M2),0,'ko','MarkerSize',7,'MarkerFaceColor','y')
plot(M1/(M1+M2),0,'ko','MarkerSize',5,'MarkerFaceColor','g')

% plot Lagrange points and labels
plot(LP(:,1),LP(:,2),'k+')
labels = {'L1', 'L2', 'L3', 'L4', 'L5'}
text(LP(:,1)-.2,LP(:,2),labels)

% plot title, limits, etc.
title(sprintf('m_1 = %.2f   m_2 = %.2f',M1, M2));
xlabel('x')
ylabel('y')
axis square
axis equal
xlim([-1.8 1.8])
ylim([-1.8 1.8])

% print plot to file
set(gcf, 'PaperPosition', [0 0 5 5]);
print ('-dpdf','Lagrange.pdf');
print ('-dpng','Lagrange.png','-r100');