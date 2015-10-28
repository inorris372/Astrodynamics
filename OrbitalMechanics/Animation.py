from matplotlib import pyplot as plt
from Rotation import rotate

__author__ = 'Ian'


############  OutputFcn2:  Animation  ###################################
#
# the output function for making the animation
def OutputFcn2(t, y, dy, flag):
    # don't stop
    stop = False
    # only after sucessfull steps
    if not flag:
        return

    else:
        NP = len(y) / 2  # number of particles (assumes 2D)

        if rotatingFrame:

            XYrot = rotate(y, t, -1)

            if leaveTrail:
                for ii in range(0, NP - 1):
                    trail(:, 2 * ii) = [XYrot(2 * ii - 1);
                    trail(1:end - 1, 2 * ii - 1)]
                    trail(:, 2 * ii + 1) =   [XYrot(2 * ii);
                    trail(1:end - 1, 2 * ii)]

                    plt.plot(trail(:, 2 * ii - 1), trail(:, 2 * ii), '-')
                    plt.hold(True)

            for ii in range (0:NP-1):
            plt.plot(XYrot(2 * ii), XYrot(2 * ii + 1), '.')
            plt.hold(True)
        plt.plot(LP(:, 1), LP(:, 2), 'k+')
        plt.plot(R10(1), R10(2), 'ro', 'MarkerFaceColor', 'r')
        plt.plot(R20(1), R20(2), 'go', 'MarkerFaceColor', 'g')
    else:
    if leaveTrail:
        for ii in range (1:NP):
        trail(:, 2 * ii - 1) = [y(2 * ii - 1);
        trail(1:end - 1, 2 * ii - 1)]
        trail(:, 2 * ii) = [y(2 * ii);
        trail(1:end - 1, 2 * ii)]

        plt.plot(trail(:, 2 * ii - 1), trail(:, 2 * ii), '-')
        plt.hold(True)


for ii in range (0, NP-1):
plt.plot(y(2 * ii), y(2 * ii + 1), '.')
plt.hold(True)
R1 = rotate(R10, t, 1)
R2 = rotate(R20, t, 1)

plt.plot(R1(1), R1(2), 'ro', 'MarkerFaceColor', 'r')
plt.plot(R2(1), R2(2), 'go', 'MarkerFaceColor', 'g')

plt.axis('equal')
plt.axis(plotLimits * [-1 1 - 1 1])
plt.hold(False)
pause(animateDelay), drawnow  # flush plotting commands
return stop
