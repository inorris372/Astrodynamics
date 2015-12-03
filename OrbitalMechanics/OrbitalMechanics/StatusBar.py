from waitbar import bar
from crtbpRKN1210 import tEnd, wait
__author__ = 'Ian'


##########  OutputFcn1:  Status Bar  ####################################
#
# the output function
def OutputFcn1(t, y, dy, flag):  # ok
    # don't stop
    stop = False
    # only after sucessfull steps
    if not flag:
        stop = True
    else:
        bar(t / tEnd, wait)
    return stop
