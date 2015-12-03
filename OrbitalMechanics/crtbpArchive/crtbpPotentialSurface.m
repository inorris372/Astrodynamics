% crtbpPotentialSurface.m
%
% Plots surface rendering of pseudo-potential for the circular, restricted 
% 3-body problem evaluated in the rotating reference frame tied to m1 & m2.
%
% Assumes G = 1 and r = 1, where r  is the distance between masses m1 & m2.
% The origin of the coordinate system is placed at the center-of-mass point
%
% Uses 'phong' lighting and raytracing to render the surface
%
% Program dependencies:
%    crtbpPotential.m - returns pseudo-potential
%
% Matlab-Monkey.com  10/10/2013



%%%%%%%%%%%%%%%%%%%%%  Parameters and Initialization  %%%%%%%%%%%%%%%%%%%%%

M1 = 1;       % mass 1
M2 = 0.1;     % mass 2
M = M1 + M2;  % total mass


P = 2*pi * sqrt(1 / M);     % period from Kepler's 3rd law
omega0 = 2*pi/P;            % angular velocity of massive bodies

[X,Y] = meshgrid(-1.5:0.05:1.5);    % grid of (x,y) coordinates
U = crtbpPotential(M1, M2, X, Y);   % calculate potential on grid points

figure

% define a custom, gray color map to render surface plot
map = ones(2 , 3)*.7;
colormap(map);



%%%%%%%%%%  Left figure will show potential from directly above

subplot(1,2,1)

surf(X,Y,U,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')

title(sprintf('m_1 = %.1f   m_2 = %.1f',M1, M2));
axis square
axis equal
axis tight
zlim([-3 -1]);
view(90,90)     % viewing angle straight overhead
camlight left



%%%%%%%%%%  Right figure will show potential tilted 30 degrees

subplot(1,2,2)

surf(X,Y,U,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')

title(sprintf('m_1 = %.1f   m_2 = %.1f',M1, M2));
axis tight
axis square
axis equal
zlim([-3 -1]);
view(90,30)     % viewing angle inclined by 30 degrees
camlight left



%%%%%%%%%%  Print to file

set(gcf, 'PaperPosition', [0 -.5 8 4]);
print ('-dpdf','potential.pdf');
print ('-dpng','potential.png','-r100');



