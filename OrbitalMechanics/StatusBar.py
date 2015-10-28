__author__ = 'Ian'


##########  OutputFcn1:  Status Bar  ####################################
#
# the output function
def OutputFcn1(t, y, dy, flag):  # ok
    # don't stop
    stop = False
    # only after sucessfull steps
    if not flag:
        return
    else:
        wait = waitbar(t / tEnd, wait)
    return stop
