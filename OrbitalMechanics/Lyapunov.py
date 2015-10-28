from crtbpRKN1210 import integrator
import matplotlib.pyplot as plt
import math
from LagrangePoints import lpoints

__author__ = 'Ian'

# Numerically integrates two test particles in the circular restricted
# three-body problem using the Runge-Kutta-Nystrom 12/10 integrator to
# calculate the Lyapunov exponent.
#
# Program dependences (required to run this code):
#   crtbpRKN1210.m - Integrates 3-body problem
#   rkn1210.m      - This code may be found on the MATLAB file and was
#                    written by Rody P.S. Oldenhuis.
#
# MATLAB-Monkey.com   10/6/2013


nPeriods = 100 # number of orbital periods to run simulation

M1 = 1         # mass 1
M2 = 0.001     # mass 2
M = M1 + M2    # total mass

P = 2*math.pi * math.sqrt(1 / M)   # period from Kepler's 3rd law
omega = 2*math.pi/P           # angular velocity of massive bodies


times = [0, nPeriods*P]   # set integration limits

R = 1                    # separation between masses must be 1
r1 = -R*M2/M             # x coordinate of M1
r2 = R*M1/M              # x coordinate of M2


LP = lpoints(M2/M)  # get positions of Lagrange points


##########  Set initial conditions to match Ceres *** No Chaos
P2 = P*4.60/11.86      # ratio of Mars's period to Jupiter's
e = 0.079              # eccentricity
a = R * (P2/P)**(2/3)  # calculate semimajor axis from period
x0 = a*(1-e)           # initial position at perihelion
y0 = 0
vx0 = 0                # initial velocity at perihelion
vy0 = math.sqrt(M1*(1+e)/x0)


##########  Set initial conditions for 2:3 resonance *** evidence of chaos
# P2 = P*2/3             # 2:3 resonance
# e = 0.56               # eccentricity
# a = R * (P2/P)**(2/3)  # calculate semimajor axis from period
# x0 = a*(1+e)           # initial position
# y0 = 0
# vx0 = 0                # initial velocity
# vy0 = math.sqrt(M1*(1-e)/x0)

##########  Set initial conditions for 7:4 resonance
# P2 = P*7/4             # 7:4 resonance
# e = 0.1                # eccentricity
# a = R * (P2/P)##(2/3)  # calculate semimajor axis from period
# x0 = a*(1-e)           # initial position
# y0 = 0
# vx0 = 0                # initial velocity
# vy0 = math.sqrt(M1*(1+e)/x0)

##########  Set initial conditions 0.025 to right of L4 point
# lp = 4
# x0 = LP(lp,1) + 0.024
# y0 = LP(lp,2)
# v = math.sqrt(x0**2+y0**2)*omega
# th = math.atan2(y0,x0)
# vx0 = v * math.cos(th+math.pi/2)
# vy0 = v * math.sin(th+math.pi/2)


NP = 2                 # number of test particles = 2


##########  Set plotting flags for integrator
rotatingFlag = True
animateFlag = False
trailFlag = True
PoincareFlag = False
flags = [rotatingFlag, animateFlag, trailFlag, PoincareFlag]  # plotting flags


##########  Integrate

vals = integrator([M1, M2], [x0, y0, x0, y0-1e-8], [vx0, vy0, vx0, vy0], times, flags)

t = vals[:,0]
pos = vals[:,1:2*NP]
vel = vals[:,2*NP+1:4*NP]



##########  Plot Lyapunov exponents as a function of integration time

plt.subplot(1,2,1)
d0 = math.sqrt((pos[0,0]-pos[0,2])**2 + (pos[0,1]-pos[0,3])**2)
d = math.sqrt((pos[:,0]-pos[:,2])**2 + (pos[:,1]-pos[:,3])**2)
gamma = math.log(d/d0)/t
plt.loglog(t,gamma,'b-')
plt.ylabel('\lambda')
plt.xlabel('time')
plt.title(print('Lyapunov Exponent   (m_2/m_1 = %.3f  P/P_0 = %.3f)',M2/M1, P2/P))

plt.subplot(1,2,2)
plt.semilogy(t,d)
plt.xlabel('t')
plt.ylabel('\Delta r')