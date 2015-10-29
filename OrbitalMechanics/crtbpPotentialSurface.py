import matplotlib as mpl
import matplotlib.cm as cm
import math
from matplotlib.pyplot import *
import numpy as np
from CrtbpPotential import potential

__author__ = 'Ian'

""" Plots surface rendering of pseudo-potential for the circular, restricted
3-body problem evaluated in the rotating reference frame tied to m1 & m2.

Assumes G = 1 and r = 1, where r  is the distance between masses m1 & m2.
The origin of the coordinate system is placed at the center-of-mass point

Uses 'phong' lighting and raytracing to render the surface

Program dependencies:
   crtbpPotential.m - returns pseudo-potential
Matlab-Monkey.com  10/10/2013
"""

##############  Parameters and Initialization  ###########

M1 = 1        # mass 1
M2 = 0.1      # mass 2
M = M1 + M2   # total mass


P = 2*math.pi*(1 / M)**.5     # period from Kepler's 3rd law
omega0 = 2*math.pi/P            # angular velocity of massive bodies

x = np.arange(-1.5, 1.5, .05)
y = np.arange(-1.5, 1.5, .05)
X,Y = np.meshgrid(x,y)    # grid of (x,y) coordinates
U = potential(M1, M2, X, Y)   # calculate potential on grid points


# define a custom, gray color map to render surface plot
map = np.ones(2, 3)*.7
fig = figure()
ax = fig.gca(projection='3d')
mycmap = cm.ma..to_rgba(map)

surf = ax.plot_surface(X, Y, U, rstride=1, cstride=1, facecolors=C, antialiased=True)
norm = mpl.colors.Normalize(vmin=-20, vmax=10)
cmap = cm.hot
x = 0.3

m = cm.ScalarMappable(norm=norm, cmap=cmap)
print m.to_rgba(x)

###########  Left figure will show potential from directly above

subplot(1,2,1)

surf(X,Y,U,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')

title('m_1 = %.1f   m_2 = %.1f', M1, M2)
axis('square')
axis('equal')
axis('tight')
#plt.zlim([-3, -1])
#plt.view(90,90)     # viewing angle straight overhead
#plt.camlight left



############  Right figure will show potential tilted 30 degrees

subplot(1,2,2)

surf(X,Y,U,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')

title('m_1 = %.1f   m_2 = %.1f',M1, M2)
axis('tight')
axis('square')
axis('equal')
#plt.zlim([-3 -1]);
#plt.view(90,30)     # viewing angle inclined by 30 degrees
#camlight left

show()
#############  Print to file

#set(gcf, 'PaperPosition', [0, -.5, 8, 4])
print ('-dpdf','potential.pdf')
print ('-dpng','potential.png','-r100')
