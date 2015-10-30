import inspect
from StatusBar import OutputFcn1 as out1
from scipy import integrate
from Derivatives import derivatives
from Rotation import rotate
import time
from StatusBar import OutputFcn1
from Animation import OutputFcn2
from LagrangePoints import lpoints
from scipy import *
from rkn1210 import rkn1210
from waitbar import bar

__author__ = 'Ian'

# Numerically integrate a test particle in the circularly restricted
# three-body problem using the Runge-Kutta-Nystrom 12/10 integrator.
#
# Massive bodies M1 and M2 follow circular orbits around their center of
# mass point, set to be the origin. The third test particle does affect the
# motion of the two massive bodies.
#
# We assume G = 1 and R = 1 (distance between the primary bodies)
#
# This function is designed to be called by another function which sets
# the initial conditions.
#
# This function calls the following functions:
#   - rkn1210.m (written by Rody P.S. Oldenhuis, MATLAB File Exchange)
#   - LagrangePoints.m
#
# Passed parameters:
#   masses - 1x2 matrix containing m1 and m2
#   pos0   - 1x2N matrix cointaining pairs of (x,y) coordinates for
#            the initial conditions of the test particles. If one test
#            particle is used, pos0 will be a 1x2 matrix containing (x0, y0)
#            For 2 test particles, posq will contain (x10, y10, x20, y20)...
#   vel0   - 1x2N matrix containing pairs of (vx, vy) values for each of
#            the initial velocities of the test particles.
#   times  - 1x2 array containing begin and end times [tBegin, tEnd]
#   flag   - 1x3 array containing plotting options (flags are optional):
#            flag(1) = true to plot and return values in rotating frame
#            flag(2) = true to animate the orbit as it integrates. Cool to
#                      watch, but it slows the calculate way down.
#            flag(3) = true to leave a trail behind the test particles when
#                      animation is turned on.
#
# Returned arrays:
#   t      - Tx1 array returning times used in integration, where T =
#            number of time steps
#   pos    - Tx2N array returning (x,y) coordinate pairs for each test
#            particle
#   vel    - Tx2N array returning (vx,vy) velocity component pairs for each
#            test particle
#
#
# MATLAB-Monkey.com   10/6/2013

# global variables used in other functions listed below

global tEnd, wait, mu1, mu2, R10, R20, omega0, plotLimits, LP, leaveTrail, rotatingFrame, animateDelay


def integrator(masses, pos0, vel0, times, flag):
    global tEnd, wait, mu1, mu2, R10, R20, omega0, plotLimits, LP, leaveTrail, rotatingFrame
    # return returnVals
    (nargin, varargs, keywords, defaults) = inspect.getargspec(integrator)
    if nargin < 4:
        print('crtbpRKN1210 requires at least 4 parameters')
        exit()



    ##########                                                       ##########
    ##########           Initialization and Parameters               ##########
    ##########                                                       ##########



    ####################  Flags Controlling Output  ###########################
    animate = False
    rotatingFrame = False

    if nargin >= 5:
        rotatingFrame = flag(1)  # if true, transform to rotating frame
        # if false, display inertial frame

        animate = flag(2)  # if true, animate motion
        # if false, just display final trajectory

        leaveTrail = flag(3)  # if true and animation is on, past positions will
        # not be erased.  If false and animation is on,
        # no 'trail' of past positions will be left

    progressBar = False  # if true, show progress bar

    animateDelay = 0.0  # delay (in seconds) between frames of animation.
    # smaller the value, faster the animation will play

    leaveTrail = flag(3)  # if true and animation is on, past positions will
    # not be erased.  If false and animation is on,
    # no 'trail' of past positions will be left

    trailLength = 700  # if leaveTrail is true, this sets the length of
    # trail

    #########################  System Parameters  ############################


    NP = len(pos0) / 2  # number of test particles.  Assumes 2D

    M1 = masses(1)  # mass 1
    M2 = masses(2)  # mass 2
    M = M1 + M2  # total mass

    G = 1  # Gravitational Constant
    R = 1  # distance between M1 and M2 set to 1


    ################  Orbital properties of two massive bodies  ##############


    mu = G * M
    mu1 = G * M1
    mu2 = G * M2

    R10 = ([-M2 / M], [0])  # initial position of M1
    R20 = ([M1 / M], [0])  # initial position of M2

    P = 2 * math.pi * math.sqrt(R ** 3 / mu)  # period from Kepler's 3rd law
    omega0 = 2 * math.pi / P  # angular velocity of massive bodies

    plotLimits = 1.5  # limit for plotting
    tEnd = times(len(times) - 1)

    LP = lpoints(M2 / M)  # calculate Lagrange points


    ###########################   Trail Properties  ###########################

    trail = ones(trailLength, 2 * NP)
    for j in range(0, 2 * NP - 1):
        trail[:, j] = pos0(j) * ones(trailLength, 1)


    ##########                                                       ##########
    ##########                      Integrate                        ##########
    ##########                                                       ##########


    # Use Runge-Kutta 45 integrator to solve the ODE

    # initialize wait variable
    wait = True
    options = integrate.ode
    if animate:
        # options = scipy.odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'outputfcn',OutputFcn2)
        options = integrate.ode(OutputFcn2).set_integrator('vode', method='bdf', order=15)
    else:  # initialize waitbar
        if progressBar:
            while (wait):
                w = bar(0, 'integrating...')
                wait = out1(0, 0, 0, w)
                # options = scipy.odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'outputfcn',OutputFcn1)
                options = integrate.ode(OutputFcn1).set_integrator('vode', method='bdf', order=15)
        else:
            # options = scipy.odeset('RelTol', 1e-12, 'AbsTol', 1e-12)
            options = integrate.ode(OutputFcn1).set_integrator('vode', method='bdf', order=15)


        # Use Runge-Kutta-Nystrom integrator to solve the ODE
    tic = time.time()
    # [t,pos,vel] = rkn1210(@derivatives, times, pos0, vel0, options)
    [t, pos, vel, te, ye] = rkn1210(derivatives, times, pos0, vel0, options)
    toc = time.time() - tic
    DT = toc

    if rotatingFrame:
        for j in range(0, len(t) - 1):
            pos[j, :] = rotate(pos[j, :].transpose(), t(j), -1)


        # returnVals = [t, pos, vel, YE]
    return [t, pos, vel]
