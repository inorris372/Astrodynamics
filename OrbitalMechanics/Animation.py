from matplotlib import pyplot as plt
from Rotation import rotate
import crtbpRKN1210 as c
import numpy as np
__author__ = 'Ian'


############  OutputFcn2:  Animation  ###################################
#
# the output function for making the animation
def OutputFcn2(t, y, dy, flag):
    # don't stop
    stop = False
    trail = []
    fig = plt.figure()
    # only after sucessfull steps
    if not flag:
        return

    else:
        NP = len(y) / 2  # number of particles (assumes 2D)

        if c.rotatingFrame:

            XYrot = rotate(y, t, -1)

            if c.leaveTrail:
                for ii in range(0, NP - 1):
                    trail[:, 2 * ii] = np.array([XYrot(2 * ii - 1)], [trail[: -2, 2 * ii - 1]])
                    trail[:, 2 * ii + 1] =   np.array([XYrot(2 * ii)], [trail[: -2, 2 * ii]])

                    plt.plot(trail[:, 2 * ii - 1], trail[:, 2 * ii], '-')
                    plt.hold(True)

            for ii in range (0,NP-1):
                plt.plot(XYrot(2 * ii), XYrot(2 * ii + 1), '.')
                plt.hold(True)

            plt.plot(c.LP[:, 0], c.LP[:, 1], 'k+')
            plt.plot(c.R10(0), c.R10(1), 'ro', 'MarkerFaceColor', 'r')
            plt.plot(c.R20(0), c.R20(1), 'go', 'MarkerFaceColor', 'g')
        else:
            if c.leaveTrail:
                for ii in range (0, NP-1):
                    trail[:, 2 * ii] = np.array([y(2 * ii)],[trail[:, 2 * ii]])
                    trail[:, 2 * ii + 1] = np.array([y(2 * ii + 1)],[trail[:, 2 * ii + 1]])

                    plt.plot(trail[:, 2 * ii], trail[:, 2 * ii + 1], '-')
                    plt.hold(True)

            for ii in range (0, NP-1):
                plt.plot(y(2 * ii), y(2 * ii + 1), '.')
                plt.hold(True)

            r1 = rotate(c.R10, t, 1)
            r2 = rotate(c.R20, t, 1)

            plt.plot(r1(1), r1(2), 'ro', 'MarkerFaceColor', 'r')
            plt.plot(r2(1), r2(2), 'go', 'MarkerFaceColor', 'g')

        plt.axis('equal')
        plt.axis(c.plotLimits * [-1,1, - 1, 1])
        plt.hold(False)
        plt.pause(c.animateDelay) # flush plotting commands
        plt.draw()
        fig.canvas.flush_events()
    return stop
