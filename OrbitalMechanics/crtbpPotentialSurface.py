import math
from matplotlib import pyplot as plt
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


# define a custom, green color map to render surface plot

map = np.ones(2, 3)*.7

fig = plt.figure()
ax = plt.axes(projection='3d')
norm = plt.Normalize()
color = plt.cm.ScalarMappable(norm(map))
surf = ax.plot_surface(X, Y, U,linewidth = 0.3, rstride=1, cstride=1, facecolors= color, antialiased=True, alpha = 0)

x = 0.3

#m = cm.ScalarMappable(norm=norm, cmap=)


###########  Left figure will show potential from directly above

plt.subplot(1,2,1)

surf(X,Y,U,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')

plt.title('m_1 = %.1f   m_2 = %.1f', M1, M2)
plt.axis('square')
plt.axis('equal')
plt.axis('tight')
ax.view(90,90) # viewing angler straight overhead

fig.savefig("potentialTopView.jpeg")
#plt.zlim([-3, -1])
#plt.camlight left
plt.show()


############  Right figure will show potential tilted 30 degrees

plt.subplot(1,2,2)

surf(X,Y,U,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')

plt.title('m_1 = %.1f   m_2 = %.1f',M1, M2)
plt.axis('tight')
plt.axis('square')
plt.axis('equal')
#plt.zlim([-3 -1]);
ax.view(90,30)     # viewing angle inclined by 30 degrees
fig.savefig("potential30DegTilt.jpeg")
#camlight left

plt.show()
#############  Print to file

#set(gcf, 'PaperPosition', [0, -.5, 8, 4])
print ('-dpdf','potential.pdf')
print ('-dpng','potential.png','-r100')
