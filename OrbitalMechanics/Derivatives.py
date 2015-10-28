import math
import numpy as np
__author__ = 'Ian'

##########  Derivatives  ##############################################
        #
        # Function defining derivatives dx/dt and dy/dt
        # uses the global parameters omega, R10, R20, mu1, mu2 but changeth them not
        def derivatives(tf,wf):
            NP = len(wf)/2      # number of particles (assumes 2D)

            rot = np.array([math.cos(omega0*tf), -math.sin(omega0*tf)],
                   [math.sin(omega0*tf),  math.cos(omega0*tf)])

            R1 = rot * R10
            R2 = rot * R20
            dvdt = []

            for ii in range (0,NP-1):
                X = wf(2*ii-1)            # x component of ii^th particle
                Y = wf(2*ii)              # y component of ii^th particle

                r1 = ((X-R1[0])**2 + (Y-R1[1])**2)**(3/2)
                r2 = ((X-R2[0])**2 + (Y-R2[1])**2)**(3/2)

                dvdt[2*ii] = - mu1 * (X-R1[0]) / r1 - mu2 * (X-R2[0]) / r2
                dvdt[2*ii+1] =  - mu1 * (Y-R1[1]) / r1 - mu2 * (Y-R2[1]) / r2
            return dvdt