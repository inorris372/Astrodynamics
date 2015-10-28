import math
import scipy
import numpy as np

__author__ = 'Ian'


def lpoints(mu):
    # For a given value of mu, this function computes the location
    # of all five libration points for the circular restricted three
    # body problem.  It returns them as equilibrium points in R^3
    # space.  Then the output is five points each with three components.
    # Each point is a column in a matrix with three rows.  The first column is
    # L1, the second L2, and so on to the last which is L5.


    # Compute the location of the libration points

    l = 1 - mu

    LP = scipy.zeros(5, 3)

    # L1
    p_L1 = {1, 2 * (mu - l), l ** 2 - 4 * l * mu + mu ** 2, 2 * mu * l * (l - mu) + mu - l,
            mu ** 2 * l ** 2 + 2 * (l ** 2 + mu ** 2), mu ** 3 - l ** 3}
    L1roots = np.roots(p_L1)
    # initialize L1 for loop
    L1 = 0
    for i in range(0, 4):
        if (L1roots[i] > -mu) & (L1roots[i] < l):
            L1 = L1roots[i]

    LP[0, 0] = L1


    # L2
    p_L2 = {1, 2 * (mu - l), l ** 2 - 4 * l * mu + mu ** 2, 2 * mu * l * (l - mu) - (mu + l),
            mu ** 2 * l ** 2 + 2 * (l ** 2 - mu ** 2), -(mu ** 3 + l ** 3)}
    L2roots = np.roots(p_L2)
    # initialize L2 for loop
    L2 = 0
    for i in range(0, 4):
        if (L2roots[i] > -mu) & (L2roots[i] > l):
            L2 = L2roots[i]

    LP[1, 0] = L2


    # L3
    p_L3 = {1, 2 * (mu - l), l ** 2 - 4 * mu * l + mu ** 2, 2 * mu * l * (l - mu) + (l + mu),
            mu ** 2 * l ** 2 + 2 * (mu ** 2 - l ** 2), l ** 3 + mu ** 3}
    L3roots = np.roots(p_L3)
    # initialize L3 for loop
    L3 = 0
    for i in range(0, 4):
        if L3roots[i] < -mu:
            L3 = L3roots[i]

    LP[2, 0] = L3


    # L4
    LP[3, 0] = 0.5 - mu
    LP[3, 1] = math.sqrt(3) / 2


    # L5
    LP[4, 0] = 0.5 - mu
    LP[4, 1] = -math.sqrt(3) / 2

    return LP
