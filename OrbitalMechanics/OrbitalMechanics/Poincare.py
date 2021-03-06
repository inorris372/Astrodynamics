import math
from LagrangePoints import lpoints
from crtbpRK45 import RK45
import matplotlib.pyplot as plt

__author__ = 'Ian'
"""  Calculates the Poincare section for the Circular, Restricted Three-
Body Problem using the Runge-Kutta 4/5 integrator.

crtbpRK45.m

Massive bodies M1 and M2 follow circular orbits around their center of
mass point, set to be the origin. The third test particle does affect the
motion of the two massive bodies.

We assume G = 1 and R = 1 (distance between the primary bodies)

This function calls the following functions:
   - rkn45.m
   - LagrangePoints.m

MATLAB-Monkey.com   10/18/2013
"""

nPeriods = 100            # number of orbital periods to run simulation

M1 = 1                   # mass 1
M2 = 0.001               # mass 2
M = M1 + M2              # total mass

P = 2*math.pi * math.sqrt(1 / M)   # period from Kepler's 3rd law
omega = 2*math.pi/P           # angular velocity of massive bodies

times = [0, nPeriods*P]   # set integration limits

R = 1                    # separation between masses must be 1
r1 = -R*M2/M             # x coordinate of M1
r2 = R*M1/M              # x coordinate of M2


LP = lpoints(M2/M) # get positions of Lagrange points


##########  Set initial conditions for 4:7 resonance
P2 = P*4/7                   # 4:7 resonance
e = 0.1                      # eccentricity
a = R * (P2/P)**(2/3)        # calculate semimajor axis from period
x0 = a*(1-e)                 # initial position
y0 = 0
vx0 = 0                      # initial velocity
vy0 = math.sqrt(M1*(1+e)/x0) # - x0*omega;


##########  Set initial conditions 0.025 to right of L4 point
# lp = 4
# x0 = LP(lp,1) + 0.024
# y0 = LP(lp,2)
# v = sqrt(x0^2+y0^2)*omega
# th = atan2(y0,x0)
# vx0 = v * cos(th+pi/2)
# vy0 = v * sin(th+pi/2)


#########  Set initial conditions to match Ceres
# P2 = P*4.60/11.86      # ratio of Mars's period to Jupiter's
# e = 0.079              # eccentricity
# a = R * (P2/P)^(2/3)   # calculate semimajor axis from period
# x0 = a*(1+e)           # initial position at aphelion
# y0 = 0
# vx0 = 0                # initial velocity at aphelion
# vy0 = math.sqrt(M1*(1-e)/x0)



##########  Set plotting flags for integrator
PoincareFlag = True
flags = PoincareFlag  # plotting flags


##########  Integrate
(t, pos, vel, YE) = RK45([M1, M2], [x0, y0], [vx0, vy0], times, flags)



##########  Plot orbit

plt.subplot(1,2,1)
plt.plot(pos[:,0],pos[:,1],'b-')
plt.axis('equal')
plt.axis('square')

##########  Plot poincare section
plt.subplot(1,2,2)
plt.plot(YE[:,0],YE[:,2],'bo','MarkerSize',2)
plt.ylabel('\odot{x}')
plt.xlabel('x')
plt.title(print('Poincare Section   (m_2/m_1 = %.3f  P/P_0 = %.3f)',M2/M1, P2/P))
