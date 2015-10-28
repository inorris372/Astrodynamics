import math
from crtbpRKN1210 import integrator
from LagrangePoints import lpoints

__author__ = 'Ian'

# Numerically integrate a test particle in the circularly restricted
# three-body problem using the Runge-Kutta-Nystrom 12/10 integrator.
# This is a front-end user interface for the RKN-12/10 integrator. Shows
# the trajectory of the test particle in an inertial or rotating frame of
# reference. Has the option to animate the dynamics as it integrates.
#
# This code requires the following functions to run:
#   - crtbpRKN1210.m
#   - rkn1210.m (written by Rody Oldenhuis, available on MATLAB File Exchange)
#   - LagrangePoints.m
#
# Syntax for calling crtRKN1210 code:
#
# Simulate a single particle:
# vals = crtbpRKN1210([M1 M2], [x0 y0], [vx0 vy0], [tBegin tEnd], flags);
#
# % Simulate two particles particle:
# vals = crtbpRKN1210([M1 M2], [x1 y1 x2 y2], [vx1 vy1 vx2 vy2], [tBegin tEnd], flags);
#
# MATLAB-Monkey.com   10/6/2013

nperiods = 20            # number of orbital periods to run simulation

M1 = 1                   # mass 1
M2 = 0.001               # mass 2
M = M1 + M2              # total mass

P = 2*math.pi * (1 / M)**.5   # period from Kepler's 3rd law
omega = 2*math.pi/P           # angular velocity of massive bodies

times = [0, nperiods*P]   # set integration limits

R = 1                   # separation between masses must be 1
r1 = -R*M2/M            # x coordinate of M1
r2 = R*M1/M             # x coordinate of M2

LP = lpoints(M2/M) # get positions of Lagrange points


"""  Set initial conditions to match Ceres
# P2 = P*4.60/11.86       # ratio of Ceres's period to Jupiter's
# e = 0.079               # eccentricity
# a = R * (P2/P)**(2/3)   # calculate semimajor axis from period
# x0 = a*(1-e)            # initial position at perihelion
# y0 = 0
# vx0 = 0                 # initial velocity at perihelion
# vy0 = math.sqrt(M1*(1+e)/x0)
"""

"""  Set initial conditions for 2:3 resonance
# P2 = P*2/3              % 2:3 resonance
# e = 0.56                % eccentricity
# a = R * (P2/P)**(2/3)   % calculate semimajor axis from period
# x0 = a*(1+e)            % initial position
# y0 = 0
# vx0 = 0                 % initial velocity
# vy0 = math.sqrt(M1*(1-e)/x0)
"""

#  Set initial conditions for 1:2 resonance
P2 = P/2               # 1:2 resonance
e = 0.8                # eccentricity
a = R * (P2/P)**(2/3)   # calculate semimajor axis from period
x0 = a*(1+e)           # initial position
y0 = 0
vx0 = 0                # initial velocity
vy0 = (M1*(1-e)/x0)**.5

"""  Set initial conditions 0.025 to right of L4 point
# lp = 4
# x0 = LP(lp,1) + 0.024
# y0 = LP(lp,2)
# v = math.sqrt(x0**2+y0**2)*omega
# th = math.atan2(y0,x0)
# vx0 = v * math.cos(th+pi/2)
# vy0 = v * math.sin(th+pi/2)
"""

NP = 1                 # number of test particles = 1


#  Set plotting flags for integrator
rotatingFlag = True
animateFlag = True
trailFlag = True
flags = (rotatingFlag, animateFlag, trailFlag)  # plotting flags

#  Integrate

vals = integrator([M1, M2], {x0, y0}, [vx0, vy0], times, flags)
# extract time, positions and velocities of test particle(s)
t = vals[:,0]
pos = vals[:,1:2*NP]
vel = vals[:,2*NP+1:4*NP]