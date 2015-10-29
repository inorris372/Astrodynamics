import math
import numpy as np
__author__ = 'Ian'

#########  Rotation  ##################################################
#
# rotate column matrix containing (x,y)
# if dir == 1, rotate forward
# if dir == -1, rotate backward
def rotate(XYin, tin, dir):

    """

    :param XYin:
    :param tin:
    :param dir:
    """
    global omega0
    NP = len(XYin)/2
    rot = np.array([math.cos(omega0*tin), -math.sin(omega0*tin),
                 math.sin(omega0*tin), math.cos(omega0*tin)])

    if dir == -1:
        rot = rot.transpose()
    XYout = []
    for ii in range (0,NP-1):
        XYout[2*ii:2*ii+1] = (rot * XYin[2*ii:2*ii+1]).transpose()
    return XYout

