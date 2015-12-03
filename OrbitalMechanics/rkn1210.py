import inspect
import numpy as np
from statsmodels.tsa.vector_ar import output
import matlab.engine


__author__ = 'Ian'


def rkn1210(funfcn, tspan, y0, yp0, options):
    # RKN1210() is a 12th/10th order numerical integrator for ordinary
    # differential equations of the special form
    #
    #   y'' = f(t, y)                     (1)
    #
    # with initial conditions
    #
    #   y (t0) = y0                       (2)
    #   y'(t0) = yp0
    #
    # This second-order differential equation is integrated with a
    # Runge-Kutta-Nystrom method, with 17 function evaluations per step. The
    # RKN-class of integrators is especially suited for this purpose, since
    # compared to a classic Runge-Kutta integration scheme the same accuracy
    # can be obtained with less function evaluations.
    #
    # This RKN12(10) method is a very high-order method, to be used in problems
    # with *extremely* stringent error tolerances. As the name implies, the
    # error should be less than O(h^13). In verious studies, it has been shown
    # that this particular integration technique is overall more efficient for
    # ODE's of the form (1) than multi-step or extrapolation methods that give
    # the same accuracy.
    #
    # RKN1210()'s behavior is very similar MATLAB's ODE-integrator suite:
    #
    # USAGE:
    # ----------------
    #
    #   [t, y, yp] = RKN1210(funfcn, tspan, y0, yp0)
    #   [t, y, yp] = RKN1210(funfcn, tspan, y0, yp0, options)
    #
    #   [t, y, yp, exitflag, output] = RKN1210(...)
    #   [t, y, yp, TE, YE, YPE, IE, exitflag, output] = RKN1210(...)
    #
    # INPUT ARGUMENTS:
    # ----------------
    #
    #  funfcn  - definition of the second-derivative function f(t, y)
    #            (See (1)). It should accept a scalar value [t] and a
    #            column vector [y] which has the same number of elements
    #            as the initial values [y0] and [dy0] provided.
    #
    #    tspan - time interval over which to integrate. It can be a
    #            two-element vector, in which case it will be interpreted
    #            as an interval. In case [tspan] has more than two
    #            elements, the integration is carried out to all times in
    #            [tspan]. Only the values for those times are then
    #            returned in [y] and [yp].
    #
    #  y0, yp0 - initial values as in (2). Both should contain the same number
    #            of elements.
    #
    #  options - options structure, created with eng.odeset(). Used options are
    #            MaxStep, InitialStep, AbsTol, Stats, Event, OutputFcn,
    #            OutputSel, and Refine. See the help for eng.odeset() for more
    #            information.
    #
    # How to use Event and/or Output functions is described in the documentation
    # on eng.odeset(). There is one difference: RKN1210() now also passes the first
    # derivative [yp] at each step as an argument:
    #
    #   status = outputFcn(t, y, yp, flag)
    #   [value, isterminal, direction] = event(t, y, yp)
    #
    # where [t] is scalar, and [y] and [yp] are column vectors, as with f(t,y).
    #
    #
    # OUTPUT ARGUMENTS:
    # ----------------
    #
    # t, y, yp - The approximate solutions for [y] and [y'] at times [t].
    #            All are concatenated row-wise, that is
    #
    #               t  = N-by-1
    #               y  = N-by-numel(y0)
    #               y' = N-by-numel(y0)
    #
    #            with N the number of sucessful steps taken during the
    #            integration.
    #
    #  exitflag - A scalar value, indicating the termination conditions
    #             of the integration:
    #
    #             -2: a non-finite function value was encountered during the
    #                 integration (INF of NaN); the integration was stopped.
    #             -1: the step size [h] fell below  the minimum acceptable
    #                 value at some time(s) [t]; results may be inaccurate.
    #              0: nothing was done; initial state.
    #             +1: sucessful integration, normal exit.
    #             +2: integration was stopped by one of the output
    #                 functions.
    #             +3: One or more events were detected, and their
    #                 corresponding [isterminal] condition also evaluated to
    #                 [true].
    #
    # TE,YE,    - These arguments are only returned when one or more event
    #  YPE,IE     functions are used. [TE] contains the times at which events
    #             were detected. [YE] and [YPE] lists the corresponding values
    #             of the solution [y] and the first derivative [yp] at these
    #             times. [IE] contains indices to the event-functions with
    #             which these events were detected. Use a smaller value for
    #             AbsTol (in [options]) to increase the accuracy of these
    #             roots when required.
    #
    #    output - structure containing additional information about the
    #             integration. It has the fields:
    #
    #             output_h              step size (sucesful steps only) at
    #                                   each time [tn]
    #             output_rejected       amount of rejected steps
    #             output_accepted       amount of accepted steps
    #             output_delta          estimate of the largest possible
    #                                   error at each time [tn]
    #             output_message        Short message describing the
    #                                   termination conditions
    #
    #             Note that these fields contain the information of ALL
    #             steps taken, even for cases where [tspan] contains
    #             more than 2 elements.
    #
    #
    # See also ODE45, ODE86, RKN86.



    # Author:
    # Name       : Rody P.S. Oldenhuis
    # E-mail     : oldenhuis@gmail.com
    # Affiliation: LuxSpace sarl.
    #
    # please report any bugs or suggestions to oldenhuis@gmail.com



    # Based on the code for ODE86 and RKN86, also available on the MATLAB
    # FileExchange.
    #
    # The construction of RKN12(10) is described in
    # High-Order Embedded Runge-Kutta-Nystrom Formulae
    # J. R. DORMAND, M. E. A. EL-MIKKAWY, AND P. J. PRINCE
    # IMA Journal of Numerical Analysis (1987) 7, 423-430
    #
    # Coefficients obtained from
    # http://www.tampa.phys.ucl.ac.uk/rmat/test/rknint.f
    # These are also available in any format on request to these authors.



    """
    ELEMENTARY EXAMPLE
    (2-body gravitational interaction / circular orbit):

    # Equations of motion for a circular orbit in 2D
    f1 = @(t, y) -y/sqrt(y'*y)^3;
    f2 = @(t, y) [y(3:4); -y(1:2)/sqrt(y(1:2).'*y(1:2))^3];

    # RKN1210 with moderate accuracy setting
    disp('RKN1210, with AbsTol = RelTol = 1e-6:')
    options = eng.odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    tic
    [t1, y1] = rkn1210(f1, [0 1000], [1; 0], [0; 1], options);
    toc
    disp('Maximum absolute error:')
    max( (y1(:,1)-cos(t1)).^2 + (y1(:,2)-sin(t1)).^2 )

    # This is how much ODE45 will have to be tuned to achieve similar accuracy
    print('\n\n')
    disp('RKN1210, with AbsTol = RelTol = 1e-6:')
    options = eng.odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    tic
    [t2, y2] = ode45(f2, [0, 1000], [1; 0; 0; 1], options);
    toc
    disp('Maximum absolute error:')
    max( (y2(:,1)-cos(t2)).^2 + (y2(:,2)-sin(t2)).^2 )
    """


    # CHANGELOG:
    #    2013/June/11
    #        - Set 'InitialStep' to previous step on recursive calls to RKN1210
    #
    #    2013/June/07
    #        - Added value of the second derivative for succesful steps to the output structure
    #        - No memory is assigned for output variables when no output arguments are requested
    #        - Removed all support for functions defined be string; only function handles
    #        - Added a few more checks on user inputs
    #        - Updated elementary example
    #        - Minor stylistic changes
    #

    def finalize():
        nonlocal tout, yout, dyout, output, index, exitflag, have_outputFcn, num_outputFcn, \
            have_events
        # first cut off redundant zeros
        varargout = []
        if produce_output:
            tout = tout[0:index,]
            yout = yout[0:index,]
            dyout = dyout[0:index,]
            output().h = output().h[0:index,]
            output().delta = output().delta[0:index,]

        # neutral flag means all OK
        if (exitflag == 0):
            exitflag = 1

        # final call to the output functions
        if have_outputFcn:
            for kk in range (0, num_outputFcn):
                try:
                    OutputFcn[kk](t, y[OutputSel], dy[OutputSel], 'done')

                except ME:
                    ME2 = ValueError('rkn1210:OutputFcn_failure_finalization\n \
                        Output function ##1d failed to evaluate during final call.', k)
                    raise ME2

        # add message to output structure
        if produce_output:
            if exitflag == 1:
                output().message = 'Integration completed successfully.'
            elif exitflag == 2:
                output().message = 'Integration terminated by one of the output functions.'
            elif exitflag == 3:
                output().message = 'Integration terminated by one of the event functions.'
            elif exitflag == - 1:
                output().message = print(['Integration terminated successfully, but the step size\n \
                fell below the minimum allowable value of #6.6e.\n \
                for one or more steps. Results may be inaccurate.'], hmin)
            elif exitflag == - 2:
                output().message = print(['Integration unsuccessful; second derivative function\n \
                returned a NaN or INF value in the last step.'])

            # handle varargout
            if have_events:
                varargout[1] = TE
                varargout[2] = YE
                varargout[3] = YPE
                varargout[4] = IE
                varargout[5] = exitflag
                varargout[6] = output

            else:
                varargout[1] = exitflag
                varargout[2] = output

        # display stats
        if Stats == 'on':
            print(1, '\n\n Number of successful steps     : #6d\n', output().accepted)
            print(1, ' Number of rejected steps       : #6d\n', output().rejected)
            print(1, ' Number of evaluations of f(t,y): #6d\n\n', output().fevals)
            print(1, '#s\n\n', output().message)

    def detect_Events(which_event, ta, fa, tb, fb):

        nonlocal output, index, TE, YE, YPE, IE
        # Detect events
        eng = matlab.engine.start_matlab('-nojvm')
        # Use false - position method(derivative - free)
        # initialize
        opts = options
        iterations = 0
        maxiterations = 1e4
        y0 = yp
        dy0 = dyp
        tt = ta

        # prune unnecessary options, and set initial step to current step
        opts = eng.odeset(opts, 'Event', [], 'OutputFcn', [], 'Stats', 'off', 'InitialStep', h)

        # start iteration
        while (min(abs(fa), abs(fb)) > abstol):

            # make step
            iterations = iterations + 1
            ttp = tt
            tt = (0.5 * fb * ta - fa * tb) / (0.5 * fb - fa)

            # safety precaution
            if (ttp == tt):
                break

            # Evaluating the event - function at this new trial location is
            # somewhat complicated.  We need to recursively call this
            # RKN1210 - routine to get appropriate values for [y] and [dy] at
            # the new time[tt] into the event function:
            [null, Zyout, Zdyout, null, Zoutput] = rkn1210(funfcn, [ttp, tt], y0, dy0, opts)

            # set new initial step to next - to - last step of previous call
            opts = eng.odeset(opts, 'InitialStep', Zoutput().h(-2))

            # save old values for next iteration
            y0 = Zyout[-1,:].transpose()
            dy0 = Zdyout[-1,:].transpose()

            # NOW evaluate event - function with these values
            fval = Event[which_event](tt, y0, dy0)

            # keep track of number of function evaluations
            if produce_output:
                output().fevals = output().fevals + Zoutput().fevals

            # compute new step
            if (fb * fval > 0):
                tb = tt
                fb = fval
            else:
                ta = tt
                fa = fval


            # check no. of iterations
            if (iterations > maxiterations):
                raise ValueError('rkn1210:rootfinder_exceeded_max_iterations',
                    'Root could not be located within #d iterations. Exiting...', maxiterations)

        # End Regula-falsi loop

        # The zero has been found!
        if produce_output:

            # insert values into proper arrays
            TE = np.array([TE], [tt])
            YPE = np.array([YPE], [dy0])
            YE = np.array([YE], [y0])
            IE = np.array([IE], [which_event])

            # if new values won't fit, first grow the arrays
            if index > np.size(yout, axis = 0):
                grow_arrays()

            # The integrand first overshoots the zero; that's how it's
            # detected.  We want the zero to be in the final arrays, but we also
            # want them in chronological order.So, move the overshoot one
            # down, and insert the zero in its place:

            index += 1
            yout[index,:] = yout[index - 1,:]
            yout[index - 1,:] = y0.transpose()
            dyout[index,:] = dyout[index - 1,:]
            dyout[index - 1,:] = dy0.transpose()
            tout[index,:] = tout[index - 1,:]
            tout[index - 1,:] = tt
            output().h[index,:] = tt - tout[index - 1,:]
            # find zeros of Event - functions
        eng.quit()

    def grow_arrays():
        nonlocal output, tout, yout, dyout
        if produce_output:

            init = 5e2
            nans = np.isnan(init, 1)
            nans_y = np.isnan(init, len(y0)-1)

            tout = np.array([tout], [nans])
            yout = np.array([yout], [nans_y])
            dyout = np.array([dyout], [nans_y])

            output().delta = np.array([output().delta], [nans])
            output().h = np.array([output().h], [nans])
            output().d2ydt2 = np.array([output().d2ydt2], [nans_y])

    eng = matlab.engine.start_matlab('-nojvm')

    # Error traps

    # args error check
    # narginchk(4,5);
    (nargin, varargs, keywords, defaults) = inspect.getargspec(rkn1210)
    # further checks on input
    assert isinstance(funfcn, 'function_handle'), 'rkn1210:funfcn_isnt_a_function'' \
    ''\nSecond derivative f(t,y) must be given as function handle.'
    assert (not tspan) & (len(tspan)-1) >= 2, 'rkn1210:tspan_empty'' \
    ''\nTime interval [tspan] must contain at least 2 values.'
    assert all(np.diff(tspan, n = 1, axis = 0) != 0), 'rkn1210:tspan_dont_span'' \
    ''\nValues in [tspan] must all span a non-zero interval.'
    assert all(np.diff(tspan, n = 1, axis = 0) > 0) or all(np.diff(tspan, n = 1, axis = 0) < 0), 'rkn1210:tspan_must_increase'' \
    ''\nThe entries in tspan must strictly increase or decrease.'
    assert len(y0) == len(yp0), 'rkn1210:initial_values_dtypegree'' \
    ''\nInitial values [y0] and [yp0] must contain the same number of elements.'
    if nargin == 5:
        assert isinstance(options,'eng.odeset()'), 'rkn1210:options_not_struct''\
        ''\nOptions must be given as a structure created with eng.odeset().'
        # Load the coefficients

    #  c
    c = np.array([0.0e0,
                  2.0e-2,
                  4.0e-2,
                  1.0e-1,
                  1.33333333333333333333333333333e-1,
                  1.6e-1,
                  5.0e-2,
                  2.0e-1,
                  2.5e-1,
                  3.33333333333333333333333333333e-1,
                  5.0e-1,
                  5.55555555555555555555555555556e-1,
                  7.5e-1,
                  8.57142857142857142857142857143e-1,
                  9.45216222272014340129957427739e-1,
                  1.0e0,
                  1.0e0])

    # matrix A is lower triangular. It's easiest to
    # load the coefficients row-by-row:
    A = np.zeros(17)

    A[1, 0:0] = 2.0e-4

    A[2, 0:1] = np.array([2.66666666666666666666666666667e-4,
                          5.33333333333333333333333333333e-4])

    A[3, 0:2] = np.array([2.91666666666666666666666666667e-3,
                          -4.16666666666666666666666666667e-3,
                          6.25e-3])

    A[4, 0:3] = np.array([1.64609053497942386831275720165e-3,
                          0.0e0,
                          5.48696844993141289437585733882e-3,
                          1.75582990397805212620027434842e-3])

    A[5, 0:4] = np.array([1.9456e-3,
                          0.0e0,
                          7.15174603174603174603174603175e-3,
                          2.91271111111111111111111111111e-3,
                          7.89942857142857142857142857143e-4])

    A[6, 0:5] = np.array([5.6640625e-4,
                          0.0e0,
                          8.80973048941798941798941798942e-4,
                          -4.36921296296296296296296296296e-4,
                          3.39006696428571428571428571429e-4,
                          -9.94646990740740740740740740741e-5])

    A[7, 0:6] = np.array([3.08333333333333333333333333333e-3,
                          0.0e0,
                          0.0e0,
                          1.77777777777777777777777777778e-3,
                          2.7e-3,
                          1.57828282828282828282828282828e-3,
                          1.08606060606060606060606060606e-2])

    A[8, 0:7] = np.array([3.65183937480112971375119150338e-3,
                          0.0e0,
                          3.96517171407234306617557289807e-3,
                          3.19725826293062822350093426091e-3,
                          8.22146730685543536968701883401e-3,
                          -1.31309269595723798362013884863e-3,
                          9.77158696806486781562609494147e-3,
                          3.75576906923283379487932641079e-3])

    A[9, 0:8] = np.array([3.70724106871850081019565530521e-3,
                          0.0e0,
                          5.08204585455528598076108163479e-3,
                          1.17470800217541204473569104943e-3,
                          -2.11476299151269914996229766362e-2,
                          6.01046369810788081222573525136e-2,
                          2.01057347685061881846748708777e-2,
                          -2.83507501229335808430366774368e-2,
                          1.48795689185819327555905582479e-2])

    A[10, 0:9] = np.array([3.51253765607334415311308293052e-2,
                           0.0e0,
                           -8.61574919513847910340576078545e-3,
                           -5.79144805100791652167632252471e-3,
                           1.94555482378261584239438810411e0,
                           -3.43512386745651359636787167574e0,
                           -1.09307011074752217583892572001e-1,
                           2.3496383118995166394320161088e0,
                           -7.56009408687022978027190729778e-1,
                           1.09528972221569264246502018618e-1])

    A[11, 0:10] = np.array([2.05277925374824966509720571672e-2,
                            0.0e0,
                            -7.28644676448017991778247943149e-3,
                            -2.11535560796184024069259562549e-3,
                            9.27580796872352224256768033235e-1,
                            -1.65228248442573667907302673325e0,
                            -2.10795630056865698191914366913e-2,
                            1.20653643262078715447708832536e0,
                            -4.13714477001066141324662463645e-1,
                            9.07987398280965375956795739516e-2,
                            5.35555260053398504916870658215e-3])

    A[12, 0:11] = np.array([-1.43240788755455150458921091632e-1,
                            0.0e0,
                            1.25287037730918172778464480231e-2,
                            6.82601916396982712868112411737e-3,
                            -4.79955539557438726550216254291e0,
                            5.69862504395194143379169794156e0,
                            7.55343036952364522249444028716e-1,
                            -1.27554878582810837175400796542e-1,
                            -1.96059260511173843289133255423e0,
                            9.18560905663526240976234285341e-1,
                            -2.38800855052844310534827013402e-1,
                            1.59110813572342155138740170963e-1])

    A[13, 0:12] = np.array([8.04501920552048948697230778134e-1,
                            0.0e0,
                            -1.66585270670112451778516268261e-2,
                            -2.1415834042629734811731437191e-2,
                            1.68272359289624658702009353564e1,
                            -1.11728353571760979267882984241e1,
                            -3.37715929722632374148856475521e0,
                            -1.52433266553608456461817682939e1,
                            1.71798357382154165620247684026e1,
                            -5.43771923982399464535413738556e0,
                            1.38786716183646557551256778839e0,
                            -5.92582773265281165347677029181e-1,
                            2.96038731712973527961592794552e-2])

    A[14, 0:13] = np.array([-9.13296766697358082096250482648e-1,
                            0.0e0,
                            2.41127257578051783924489946102e-3,
                            1.76581226938617419820698839226e-2,
                            -1.48516497797203838246128557088e1,
                            2.15897086700457560030782161561e0,
                            3.99791558311787990115282754337e0,
                            2.84341518002322318984542514988e1,
                            -2.52593643549415984378843352235e1,
                            7.7338785423622373655340014114e0,
                            -1.8913028948478674610382580129e0,
                            1.00148450702247178036685959248e0,
                            4.64119959910905190510518247052e-3,
                            1.12187550221489570339750499063e-2])

    A[15, 0:14] = np.array([-2.75196297205593938206065227039e-1,
                            0.0e0,
                            3.66118887791549201342293285553e-2,
                            9.7895196882315626246509967162e-3,
                            -1.2293062345886210304214726509e1,
                            1.42072264539379026942929665966e1,
                            1.58664769067895368322481964272e0,
                            2.45777353275959454390324346975e0,
                            -8.93519369440327190552259086374e0,
                            4.37367273161340694839327077512e0,
                            -1.83471817654494916304344410264e0,
                            1.15920852890614912078083198373e0,
                            -1.72902531653839221518003422953e-2,
                            1.93259779044607666727649875324e-2,
                            5.20444293755499311184926401526e-3])

    A[16, 0:15] = np.array([1.30763918474040575879994562983e0,
                            0.0e0,
                            1.73641091897458418670879991296e-2,
                            -1.8544456454265795024362115588e-2,
                            1.48115220328677268968478356223e1,
                            9.38317630848247090787922177126e0,
                            -5.2284261999445422541474024553e0,
                            -4.89512805258476508040093482743e1,
                            3.82970960343379225625583875836e1,
                            -1.05873813369759797091619037505e1,
                            2.43323043762262763585119618787e0,
                            -1.04534060425754442848652456513e0,
                            7.17732095086725945198184857508e-2,
                            2.16221097080827826905505320027e-3,
                            7.00959575960251423699282781988e-3,
                            0.0e0])

    # this facilitates and speeds up implementation
    A = A.transpose()
    # Bhat (high-order b)
    Bhat = np.array([1.21278685171854149768890395495e-2,
                     0.0e0,
                     0.0e0,
                     0.0e0,
                     0.0e0,
                     0.0e0,
                     8.62974625156887444363792274411e-2,
                     2.52546958118714719432343449316e-1,
                     -1.97418679932682303358307954886e-1,
                     2.03186919078972590809261561009e-1,
                     -2.07758080777149166121933554691e-2,
                     1.09678048745020136250111237823e-1,
                     3.80651325264665057344878719105e-2,
                     1.16340688043242296440927709215e-2,
                     4.65802970402487868693615238455e-3,
                     0.0e0,
                     0.0e0])
    # BprimeHat (high-order b-prime)
    Bphat = np.array([1.21278685171854149768890395495e-2,
                      0.0e0,
                      0.0e0,
                      0.0e0,
                      0.0e0,
                      0.0e0,
                      9.08394342270407836172412920433e-2,
                      3.15683697648393399290429311645e-1,
                      -2.63224906576909737811077273181e-1,
                      3.04780378618458886213892341513e-1,
                      -4.15516161554298332243867109382e-2,
                      2.46775609676295306562750285101e-1,
                      1.52260530105866022937951487642e-1,
                      8.14384816302696075086493964505e-2,
                      8.50257119389081128008018326881e-2,
                      -9.15518963007796287314100251351e-3,
                      2.5e-2])
    # B (low-order b)
    B = np.array([1.70087019070069917527544646189e-2,
                  0.0e0,
                  0.0e0,
                  0.0e0,
                  0.0e0,
                  0.0e0,
                  7.22593359308314069488600038463e-2,
                  3.72026177326753045388210502067e-1,
                  -4.01821145009303521439340233863e-1,
                  3.35455068301351666696584034896e-1,
                  -1.31306501075331808430281840783e-1,
                  1.89431906616048652722659836455e-1,
                  2.68408020400290479053691655806e-2,
                  1.63056656059179238935180933102e-2,
                  3.79998835669659456166597387323e-3,
                  0.0e0,
                  0.0e0])
    # Bprime (low-order bprime)
    Bp = np.array([1.70087019070069917527544646189e-2,
                   0.0e0,
                   0.0e0,
                   0.0e0,
                   0.0e0,
                   0.0e0,
                   7.60624588745593757356421093119e-2,
                   4.65032721658441306735263127583e-1,
                   -5.35761526679071361919120311817e-1,
                   5.03182602452027500044876052344e-1,
                   -2.62613002150663616860563681567e-1,
                   4.26221789886109468625984632024e-1,
                   1.07363208160116191621476662322e-1,
                   1.14139659241425467254626653171e-1,
                   6.93633866500486770090602920091e-2,
                   2.0e-2,
                   0.0e0])

    #
    #
    # initialize all variables
    nargout = 0
    exitflag = 0
    pow = 1 / 12
    t0 = tspan[0]
    tfinal = tspan[len(tspan) - 1]
    t = t0
    y = y0.flatten(1)
    dy = yp0.flatten(1)
    tout = t0
    yout = y0[:].transpose()
    dyout = yp0[:].transpose()
    hmin = abs(tfinal - t) / 1e12
    f = y0[:] * np.zeros(1, 17)
    have_events = False
    halt = False
    have_outputFcn = False
    produce_output = (nargout != 0)


    # initialize output-structure
    def output():
        h = []
        fevals = 0
        rejected = 0
        accepted = 0
        delta = []
        message = 'Integration not started.'
        d2ydt2 = []

    # times might be given in reverse
    direction = 1 - 2 * ((len(tspan)-1) == 2 & tfinal < t0)

    # parse options
    if (nargin < 5):
        options = eng.odeset
        # initial options structure

    Stats = eng.odeget(options, 'Stats', 'off')
    # display statistics at the end?
    abstol = eng.odeget(options, 'AbsTol', 1e-14)  # Absolute tolerance
    reltol = eng.odeget(options, 'RelTol', 1e-7)  # Relative tolerance
    hmax = eng.odeget(options, 'MaxStep', tfinal - t)  # defaults to entire interval
    initialstep = eng.odeget(options, 'InitialStep')  # default determined later
    Event = eng.odeget(options, 'Event', [])  # defaults to no eventfunctions
    OutputFcn = eng.odeget(options, 'OutputFcn', [])  # defaults to no output functions
    num_events = 0  # initialize
    num_outputFcn = 0 # initialize
    TE = []
    TEI = []
    YE = []
    YEI = []
    YPE = []
    DYEI = []
    IE = []
    IEI = []
    # FUTURE WORK
    # NormControl = eng.odeget(options, 'NormControl', 'off')  # relax error control
    # MassMatrix = eng.odeget(options, 'Mass', [])   # Mass matrix, for problems of the...
    # ( and friends) form  M(t, y) * y = F(t, y)
    ME = BaseException
    # in case of Event / output functions, define a few extra parameters
    if (Event):

        # so we DO have to evaluate event - functions
        have_events = True
        # cast into cell - array if only one function is given(easier later on)
        if isinstance(Event, 'function_handle'):
            Event = np.array([Event])
        # number of functions provided
        num_events = len(Event)

        # Check if all are indeed function handles
        for ii in range(0, num_events):
            if not (Event[ii]) or not isinstance(Event[ii], 'function_handle'):
                raise ValueError('rkn1210:event_not_function_handles''\
                ''\nUnsupported class for event function received; event #d is of class ''#s''.\n''\
                ''RKN1210 only supports function handles.', ii, type(Event[ii]))

        # initialize TE (event times), YE (event solutions) YPE (event
        # derivs) and IE(indices to corresponding event function).Check
        # user - provided event functions at the same time

        previous_event_values = np.zeros(num_events, 1)
        for k in range (0, num_events):
            try:
                previous_event_values[k] = Event[k](t, y, dy)
            except ME:
                ME2 = RuntimeError('rkn1210:eventFcn_dont_evaluate',
                    print('Event function ##1d failed to evaluate on initial call.', k))
                raise ME2

        if produce_output:
            TE = []
            YE = []
            YPE = []
            IE = []

    Refine = 1
    if OutputFcn:

        # so we DO have to evaluate output - functions
        have_outputFcn = True
        # cast into cell - array if only one function is given(easier later on)
        if isinstance(OutputFcn, 'function_handle'):
            OutputFcn = np.array([OutputFcn])
        # number of functions provided
        num_outputFcn = len(OutputFcn)

        # Check if all are indeed function handles
        for ii in range(0, num_outputFcn):
            if not OutputFcn[ii] or not isinstance(OutputFcn[ii], 'function_handle'):
                raise ValueError('rkn1210:output_not_function_handles'' \
                    ''\nUnsupported class for output function received; output function #d is of class ''#s''.\n'' \
                    ''RKN1210 only supports function handles.', ii, type(OutputFcn[ii]))

        # WHICH elements should be passed to the output function?
        OutputSel = eng.odeget(options, 'OutputSel', [0,len(y0)-1])
        # adjust the number of points passed to the output functions by this factor
        Refine = eng.odeget(options, 'Refine', 1)

        # call each output function with 'init' flag.Also check whether
        # the user-provided output function evaluates
        for k in range (0,num_outputFcn):
            try:
                OutputFcn[k](t, y(OutputSel), dy(OutputSel), 'init')

            except ME:
                ME2 = RuntimeError('rkn1210:OutputFcn_doesnt_evaluate',
                    print('Output function ##d failed to evaluate on initial call.', k))
                raise ME2

    # Different use case: numel(tspan) > 2

    # do the calculation recursively if [tspan] has
    # more than two elements
    if len(tspan)-1 > 2:

        # the output times are already known
        tout = tspan[:]

        # call this function as many times as there are times in [tspan]
        for ii in range (0,len(tspan) - 1):

            # new initial values
            tspanI = tspan[ii:ii + 1]
            y0 = yout[len(yout) - 1,:]
            yp0 = dyout[len(dyout) - 1,:]

            # call the integrator
            if have_events:
                [toutI, youtI, dyoutI, TEI, YEI, DYEI, IEI, exitflag, outputI] = rkn1210(funfcn, tspanI, y0, yp0, options)
            else:
                [toutI, youtI, dyoutI, exitflag, outputI] = rkn1210(funfcn, tspanI, y0, yp0, options)


            # new initial step is old next - to - last step
            options = eng.odeset(options, 'InitialStep', outputI().h[-2])

            # append the solutions
            yout = np.array([yout], [youtI[-1,:]])  # ok
            dyout = np.array([dyout], [dyoutI[-1,:]])  # ok
            if have_events:
                TE = np.array([TE], [TEI])  # ok
                YEI = np.array([YE], [YEI])  # ok
                DYEI = np.array([YPE], [DYEI])  # ok
                IEI = np.array([IE], [IEI])  # ok


            # process the output
            output().h = np.array([output().h], [outputI().h])
            output().fevals = output().fevals + outputI().fevals
            output().rejected = output().rejected + outputI().rejected
            output().accepted = output().accepted + outputI().accepted
            output().delta = np.array([output().delta], [outputI().delta])

            # evaluate any output functions at each[t] in [tspan]
            if have_outputFcn:
                # evaluate all functions
                for k in range (0, num_outputFcn):
                    try:
                        halt = OutputFcn[k](toutI, youtI(OutputSel), dyoutI(OutputSel), [])

                    except ME:
                        ME2 = Exception('rkn1210:OutputFcn_failure_integration',
                            print('Output function ##1d failed to evaluate during integration.', k))
                        raise ME2


                # halt integration when requested
                if halt:
                    exitflag = 2
                    finalize()
                    return

            # should we quit?
            if exitflag == -2 or exitflag == 3:
                break
        # we're done.
        index = yout.shape[0]
        finalize()
        return

# if, for more than 2 elements in tspan

# # Initial step

    try:
        f[:, 1] = funfcn[t, y]

    except ME:
        ME2 = Exception('rkn1210:incorrect_funfcnoutput',
            print('Derivative function should return a #3.0f-element column vector.', len(y0)))
        raise ME2



    output().fevals += 1  # don't forget this one :)

    if not initialstep:
        # default initial step
        h = abstol ** pow / max(max(abs([dy.transpose(), f[:, 0].transpose()])), 1e-4)
        h = min(hmax, max(h, hmin))

    else:
        # user provided initial step
        h = initialstep

    h = direction * h  # take care of direction
    h2 = h * h   # pre - compute the square

    # # The main loop

    # Initialize all variables
    index = 1
    grow_arrays()

    # the main loop
    while (abs(t - tfinal) > 0):

        # take care of final step
        if ((t + h) > tfinal):
            h = (tfinal - t)
            h2 = h * h  # also pre - compute the square

        # Compute the second - derivative
        # NOTE: 'Vectorized' in eng.odeset() has no use; we need the UPDATED
        # function values to calculate the NEW ones, i.e., the function
        # evaluations are not independent.
        for jj in range (0,16):
            f[:, jj] = funfcn[t + c[jj] * h, y + c[jj] * h * dy + h2 * f * A[:, jj]]
            output().fevals += 1

        # check for inf or NaN
        if not any(np.isfinite(f[:])):
            exitflag = -2
        # use warning(not error) to preserve output thus far
            ValueError('rkn1210:nonfinite_values',
                ['INF or NAN value encountered during the integration.\n'' \
                ''Terminating integration...'])
            finalize()
            return
        # non - finite values

        # pre - compute the sums of the products with the coefficients
        fBphat = f * Bphat
        fBhat = f * Bhat

        # Estimate the error and the acceptable error
        delta1 = max(abs(h2 * (fBhat - f * B)))  # error ~ | Y - y |
        delta2 = max(abs(h * (fBphat - f * Bp)))  # error ~ | dot{Y} - dot{y} |
        delta = max(delta1, delta2)  # worst case error

        # update the solution only if the error is acceptable
        if (delta <= abstol) and (delta / max(np.linalg.norm(y + h * dy + h2 * fBhat), np.linalg.norm(dy + h * fBphat)) <= reltol):

            # update the new solution
            index = index + 1
            tp = t
            t = t + h
            yp = y
            y = y + h * dy + h2 * fBhat
            dyp = dy
            dy = dy + h * fBphat

            # if new values won't fit, first grow the arrays
            # NOTE: This construction is WAY better than growing the arrays on
            # EACH iteration; especially for "cheap" integrands, this
            # construction causes a lot less overhead.
            if index > yout.shape[0]:
                grow_arrays()

            # insert updated values
            if produce_output:
                tout[index,:] = t
                yout[index,:] = y.transpose()
                dyout[index,:] = dy.transpose()

                output().d2ydt2[index,:] = f[:, 1].transpose()
                output().h[index,:] = h
                output().delta[index,:] = delta
                output().accepted = output().accepted + 1


            # evaluate event - funtions
            if have_events:
                # evaluate all event - functions
                for k in range (0,num_events):
                    # evaluate event function, and check if any have changed sign
                    #
                    # NOTE: although not really necessary(event functions have been
                    # checked upon initialization), use TRY - CATCH block to produce
                    # more useful errors in case something does go wrong.

                    terminate = False
                    try:
                        # evaluate function
                        [value, isterminal, zerodirection] = Event[k](t, y, dy)

                        # look for sign change
                        if (previous_event_values[k] * value < 0):
                            # ZERODIRECTION:
                            # 0: detect all zeros(default
                            # +1: detect only INcreasing zeros
                            # -1: detect only DEcreasing zeros
                            if (zerodirection == 0) or (np.sign(value) == np.sign(zerodirection)):
                                # terminate?
                                terminate = terminate or isterminal
                                # detect the precise location of the zero
                                # NOTE: try-catch is necessary to prevent things like
                                # discontinuous event - functions from resulting in
                                # unintelligible error messages
                                try:
                                    detect_Events(k, tp, previous_event_values[k], t, value)

                                except ME:
                                    ME2 = ValueError('rkn1210:eventFcn_failure_zero',
                                        print('Failed to locate a zero for event function ##1d.', k))
                                    raise ME2

                        # save new value
                        previous_event_values[k] = value

                    except ME:
                        ME2 = ValueError('rkn1210:eventFcn_failure_integration',
                            print('Event function ##1d failed to evaluate during integration.', k))
                        raise ME2


                    # do we need to terminate?
                    if terminate:
                        exitflag = 3
                        finalize()
                        return

                    # Event functions

            # evaluate output functions
            if have_outputFcn:
                # evaluate all event - functions
                for k in range (0, num_outputFcn):
                    # evaluate output function
                    try:
                        halt = OutputFcn[k](t, y[OutputSel], dy[OutputSel], [])

                    except ME:
                        ME2 = ValueError('rkn1210:OutputFcn_failure_integration',
                            print('Output function ##1d failed to evaluate during integration.', k))
                        raise ME2
                    # halt integration when requested
                    if halt:
                        exitflag = 2
                        finalize()
                        return
            # Output functions

        # rejected step: just increase its counter
        else:
            output().rejected += 1

        # accept or reject step
        # adjust the step size
        if (delta != 0):
            h = np.sign(h) * min(hmax, 0.9 * abs(h) * (abstol / delta) ** pow)
            # Use[Refine] - option when output functions are present
            if have_outputFcn:
                h = h / Refine
            # pre - compute the square
            h2 = h * h
            # check the new stepsize
            if (abs(h) < hmin):
                exitflag = -1
                # use warning to preserve results thus far
                ValueError('rkn1210:stepsize_too_small',
                ['Failure at time t = #6.6e: \n'' \
                ''Step size fell below the minimum acceptable value of #6.6d.\n'' \
                ''A singularity is likely.'], t, hmin)
        
            # adjust step - size

    # main loop
    eng.quit()
    # if the algorithm ends up here, all was ok
    finalize()





