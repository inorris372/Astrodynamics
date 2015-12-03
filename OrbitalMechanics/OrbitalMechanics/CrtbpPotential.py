__author__ = 'Ian'

""" Pseudo-potential for the circular, restricted 3-body problem evaluated in
the rotating reference frame tied to m1 & m2
Assumes G = 1 and r = 1, where r  is the distance between masses m1 & m2
The origin of the coordinate system is placed at the center of mass point

passed parameters:
   m1 = mass of object 1 (typically set m1 + m2 = 1)
   m2 = mass of object 2
   (x, y) = coordinates of position to evaluate potential. x and y may be
      single values, or arrays.

returned value = pseudo-potential at position (x,y)
  if x & y are arrays, crtbpPotential will return an array

Matlab-Monkey.com  10/10/2013
"""


def potential (m1, m2, x, y):
    omega = (m1 + m2) ** .5  # angular velocity of massive bodies

    x1 = -m2 / (m1 + m2)  # x coordinate of m1
    x2 = 1 + x1  # x coordinate of m2

    r1 = ((x - x1) ** 2 + y ** 2) ** .5  # distance from m1 to (x,y)
    r2 = ((x - x2) ** 2 + y ** 2) ** .5  # distance from m2 to (x,y)

    rVal = -omega ** 2 / 2 * (x ** 2 + y ** 2) - m1 / r1 - m2 / r2
    return rVal
