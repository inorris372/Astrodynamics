from waitbar import bar
__author__ = 'Ian'


##########  OutputFcn1:  Status Bar  ####################################
#
# the output function

def OutputFcn1(t, y, dy, flag):  # ok
    global tEnd, wait
    # don't stop
    stop = False
    # only after sucessfull steps
    if not flag:
        stop = True
    else:
        bar(t / tEnd, wait)
    return stop
