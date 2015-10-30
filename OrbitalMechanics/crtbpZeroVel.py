import math
import numpy as np
from CrtbpPotential import potential
import matplotlib.pyplot as plt
from LagrangePoints import lpoints

__author__ = 'Ian'
G = 1  # Gravitational Constant
M1 = 1  # mass 1
M2 = 0.1  # mass 2
M = M1 + M2  # total mass

R = 1  # distance between M1 and M2 set to 1


#############  Orbital properties of two massive bodies  #############
mu = G * M
mu1 = G * M1
mu2 = G * M2

R10 = ([-M2 / M], [0])  # initial position of M1
R20 = ([M1 / M], [0])  # initial position of M2

P = 2 * math.pi * math.sqrt(R ** 3 / mu)  # period from Kepler's 3rd law
omega0 = 2 * math.pi / P  # angular velocity of massive bodies

(X, Y) = np.meshgrid(np.linspace(-2, 2, 100))
U = potential(mu1, mu2, X, Y)

U0 = potential(mu1, mu2, 0, 1)

m = 2
#map = scipy.ones(m, 3) * 8



level = [-2 - 1.96 - 1.9 - 1.7 - 1.68 - 1.64]

for j in range(0, 5):
    plt.subplot(2, 3, j)

    #  contourf(X,Y,U,[-5:.5:-1])
    # added to level[j] in original code [0, -0.0001]
    plt.contourf(X, Y, U, level[j])
    plt.hold(True)
    plt.plot(-M2 / (M1 + M2), 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'y')
    plt.plot(M1 / (M1 + M2), 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'g')
    plt.axis([min(X), max(X), min(Y), max(Y)])
    plt.title(print('C_J = %.2f', -2 * level[j]))
    plt.axis('equal')
    plt.axis('square')
    plt.hold(True)

    LP = lpoints(M2 / (M1 + M2))
    plt.plot(LP[:, 1], LP[:, 1], 'k+')

# set(gcf, 'PaperPosition', [0,-1,8,4])
print('-dpng', 'zeroVelocity.png', '-r100')
