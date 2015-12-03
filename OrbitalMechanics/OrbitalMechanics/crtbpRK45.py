import inspect
import math
import time
from Derivatives import derivatives
from waitbar import bar
from StatusBar import OutputFcn1 as out1
import scipy.integrate as i
__author__ = 'Ian'

""" Numerically integrate a test particle in the circularly restricted
three-body problem using the Runge-Kutta 4/5 integrator.

Massive bodies M1 and M2 follow circular orbits around their center of
mass point, set to be the origin. The third test particle does affect the
motion of the two massive bodies.

To run this code, rkn1210.m must be downloaded from the MATLAB file
Exchange.  Thanks to Rody P.S. Oldenhuis for his File Exchange code.

We assume G = 1 and R = 1 (distance between the primary bodies)

This function is designed to be called by another function which sets
the initial conditions.

This function calls the following functions:
   - rkn45.m
   - LagrangePoints.m

Passed parameters:
   masses - 1x2 matrix containing m1 and m2
   pos0   - 1x2N matrix cointaining pairs of (x,y) coordinates for
            the initial conditions of the test particles. If one test
            particle is used, pos0 will be a 1x2 matrix containing (x0, y0)
            For 2 test particles, posq will contain (x10, y10, x20, y20)...
   vel0   - 1x2N matrix containing pairs of (vx, vy) values for each of
            the initial velocities of the test particles.
   times  - 1x2 array containing begin and end times [tBegin, tEnd]
   flag   - 1x1 array containing plotting options (flags are optional):
            flag(1) = true to create Poincare section

Returned arrays:
   t      - Tx1 array returning times used in integration, where T =
            number of time steps
   pos    - Tx2N array returning (x,y) coordinate pairs for each test
            particle
   vel    - Tx2N array returning (vx,vy) velocity component pairs for each
            test particle
   YE     - matrix containing positions and velocities for Poincare
            Section

 MATLAB-Monkey.com   10/18/2013
"""

def RK45(masses, pos0, vel0i, times, flag):

    # ---------  Events - Poincare Section  --------- #

    def events(tin,y,dy):
        # Locate the time when y1 passes through zero in an
        # increasing direction

        value = y(2)     # Detect when particle passes y1 = 0
        isterminal = 0   # Don't stop the integration
        direction = -1   # Positive direction only
        return (value, isterminal, direction)


    (nargin, varargs, keywords, defaults) = inspect.getargspec(RK45)

    if nargin < 4:
        print('crtbpRK45 requires at least 4 parameters')
        exit()



##########                                                       %%%%%%%%%%
##########           Initialization and Parameters               %%%%%%%%%%
##########                                                       %%%%%%%%%%

###################  Flags Controlling Output  ###################

    if nargin >= 5:
        Poincare = flag[0]      # if true, create Poincare section
    else:
        Poincare = False

###################  System Parameters  #############

    NP = len(pos0)/2  # number of test particles.  Assumes 2D

    M1 = masses(1)   # mass 1
    M2 = masses(2)   # mass 2
    M = M1 + M2      # total mass

    G = 1            # Gravitational Constant
    R = 1            # distance between M1 and M2 set to 1

###################  Orbital properties of two massive bodies  ############

    mu = G*M
    mu1 = G*M1
    mu2 = G*M2

    R10 = [-M2/M, 0]            # initial position of M1
    R20 = [M1/M, 0]             # initial position of M2

    P = 2*math.pi * math.sqrt(R^3 / mu)   # period from Kepler's 3rd law
    omega0 = 2*math.pi/P             # angular velocity of massive bodies


# calculate velocity components in rotating frame
    vel0 = []
    vm = math.sqrt(vel0i[0]**2 + vel0i[1]**2)
    vxu = vel0i[0]/vm
    vyu = vel0i[1]/vm
    vr = omega0 * math.sqrt(pos0[0]**2+pos0[1]**2)
    vnew = vm-vr
    vel0[0] = vxu * vnew
    vel0[1] = vyu * vnew





# --------                                                       -------- #
# --------                      Integrate                        -------- #
# --------                                                       -------- #


# Use Runge-Kutta 45 integrator to solve the ODE


    if Poincare:
        options = i.ode(out1).set_integrator('vode', method='bdf', order=15)

        # options = .odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events',@events,'outputfcn',@OutputFcn1)
    else:
        options = i.ode(out1).set_integrator('vode', method='bdf', order=15)

        # options = sp.integrate.odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'outputfcn',@OutputFcn1)


# initialize waitbar
    wait = True
    while (wait):
                w = bar(0, 'integrating...')
                wait = out1(0, 0, 0, w)


# Use Runge-Kutta-Nystrom integrator to solve the ODE
    tic = time.time()
    [t, q, TE, YE] = i.odeint(derivatives, pos0, times[0],
                           #ode45(@derivatives, times, [pos0 vel0]', options)
    toc = time.time() - tic
    DT = toc

    pos = q[:, 0:1]
    vel = q[:, 2:3]
    return (t, pos, vel, YE)