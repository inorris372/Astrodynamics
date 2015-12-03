from LagrangePoints import lpoints
from CrtbpPotential import potential
from crtbpZeroVel import X, Y, U, mu1, mu2
import matplotlib.pyplot as plt

__author__ = 'Ian'

# Plots Lagrange points and zero-velocity curves that pass through them
# for the circular, restricted 3-body problem.
#
# Assumes G = 1 and r = 1, where r  is the distance between masses m1 & m2.
# The origin of the coordinate system is placed at the center-of-mass point
#
# Program dependencies:
#    crtbpPotential.m - returns pseudo-potential
#    lagrangePoints.m - returns 5x3 array containing (x,y,z) coordinates
#                       of the Lagrange points for given values of m1, m2
#
# Matlab-Monkey.com  10/10/2013

# --------------  Parameters and Initialization  ----------------- #

M1 = 1       # mass 1
M2 = 0.1     # mass 2
M = M1 + M2  # total mass

# P = 2*math.pi * math.sqrt(R**3 / mu)   # period from Kepler's 3rd law
# omega0 = 2*math.pi/P            # angular velocity of massive bodies

# find Lagrange points and the pseudo-potential
LP = lpoints(M2/(M1+M2))
LP1_level = potential(mu1, mu2, LP[0, 0], LP[0, 1])
LP2_level = potential(mu1, mu2, LP[1, 0], LP[1, 1])
LP3_level = potential(mu1, mu2, LP[2, 0], LP[2, 1])


# -------------------     Plotting     ----------------------- #

# plot zero-velocity curves that run through L1, L2, L3
plt.contour(X, Y, U, [LP1_level, LP2_level, LP3_level])
plt.hold(True)

# plot m1 (yellow dot) and m2 (green dot)
plt.plot(-M2/(M1+M2), 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'y')
plt.plot(M1/(M1+M2), 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g')

# plot Lagrange points and labels
plt.plot(LP[:, 0], LP[:, 1], 'k+')
labels = {'L1', 'L2', 'L3', 'L4', 'L5'}
plt.text(LP[:, 0]-.2, LP[:, 1], labels)

# plot title, limits, etc.
plt.title(print('m_1 = %.2f   m_2 = %.2f', M1, M2))
plt.xlabel('x')
plt.ylabel('y')
plt.axis('square')
plt.axis('equal')
plt.axis([-1.8, 1.8, -1.8, 1.8])

# print plot to file
# set(gcf, 'PaperPosition', [0, 0, 5, 5])
print('-dpdf', 'Lagrange.pdf')
print('-dpng', 'Lagrange.png', '-r100')
