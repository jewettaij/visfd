#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import special


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """
    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)



def main(argv):

    A = 2.7
    B = 1.7
    a = 4.0
    b = 5.0

    try:
        i=1
        if sys.argv[1] == '-gauss':
            if i+2 >= len(sys.argv):
                raise InputError('Error: Expected 2 numeric arguments following \"'+sys.argv[1]+'\"\n')
            A = float(sys.argv[i+1])
            a = float(sys.argv[i+2])
            m = 2.0
            n = 2.0
            B = 0.0
            b = 1.0
            width = a
            r_i = np.arange(-ceil(4.0*width), ceil(4.0*width), 1.0)
            # Old method: Using a Gaussian kernel
            #p_i = np.array([A*exp(-0.5*(r/a)**2) for r in r_i])
            # New method: Using a discrete Gaussian kernel
            # Unfortunatley, descrete Gaussians are numerically unstable
            # for large arguments, so we have to loop through the entries
            # and check whether the arguments for special.iv() are too large.
            p_i = np.zeros(r_i.size)
            for i in range(0, r_i.size):
                r = r_i[i]
                if a == 0:
                    if r == 0:
                        p_i[i] = 1.0
                    else:
                        p_i[i] = 0.0
                elif ((a <= 10.0) and (abs(r) <= 10)):
                    p_i[i] = exp(-a*a)*special.iv(abs(r), a*a)
                else:
                    # otherwise, specia.iv() is probably numerically unstable
                    # instead compute this the ordinary way
                    p_i[i] = exp(-0.5*(r/a)**2) / sqrt(2*pi*a*a)
                if r == 0:
                    p_0 = p_i[i]
            # Now take into account the "A" coefficient (for normalization).
            p_i *= A / p_0
            #plt.plot(r_i, p_i)
            plt.step(r_i+0.5, p_i)
            plt.show()

        elif sys.argv[1] == '-ggauss':
            if i+4 >= len(sys.argv):
                raise InputError('Error: Expected 4 numeric arguments following \"'+sys.argv[1]+'\"\n')
            A = float(sys.argv[i+1])
            a = float(sys.argv[i+2])
            #B = 0.0
            #b = 1.0
            # if the user manually specified (m and n)
            # then we assume the user wants us to plot this instead:
            #   h(r) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)
            m = float(sys.argv[i+3])
            #n = float(sys.argv[i+4])
            width = a
            r_i = np.arange(-ceil(4.0*width), ceil(4.0*width), 1.0)
            p_i = np.array([A*exp(-((abs(r)/a)**m)) for r in r_i])
            #plt.plot(r_i, p_i)
            plt.step(r_i+0.5, p_i)
            plt.show()

        elif sys.argv[1] == '-dog':
            if i+4 >= len(sys.argv):
                raise InputError('Error: Expected 4 numeric arguments following \"'+sys.argv[1]+'\"\n')
            A = float(sys.argv[i+1])
            B = float(sys.argv[i+2])
            a = float(sys.argv[i+3])
            b = float(sys.argv[i+4])
            width = max(a,b)
            r_i = np.arange(-ceil(4.0*width), ceil(4.0*width), 1.0)
            # Old method: Using a Gaussian kernel
            #p_i = np.array([A*exp(-0.5*(r/a)**2) -
            #                B*exp(-0.5*(r/b)**2) for r in r_i])
            # New method: Using a discrete Gaussian kernel.
            # I will use a for-loop here because it's really more messy to try
            # and do this with a python-style list-comprehension one-liner.
            # Either way, this is ugly code, and I'm too lazy to clean it.

            # Calculate the first Guassian kernel
            pa_i = np.zeros(r_i.size)
            for i in range(0, r_i.size):
                r = r_i[i]
                if a == 0:
                    if r == 0:
                        pa_i[i] = 1.0
                    else:
                        pa_i[i] = 0.0
                elif ((a <= 10.0) and (abs(i) <= 10)):
                    pa_i[i] = exp(-a*a)*special.iv(abs(r), a*a)
                else:
                    # otherwise, specia.iv() is probably numerically unstable
                    # instead compute this the ordinary way
                    pa_i[i] = exp(-0.5*(r/a)**2) / sqrt(2*pi*a*a)
                if r == 0:
                    pa_0 = pa_i[i]
            # Now take into account the "A" coefficient (for normalization).
            pa_i *= A / pa_0

            # Calculate the second Guassian kernel
            pb_i = np.zeros(r_i.size)
            for i in range(0, r_i.size):
                r = r_i[i]
                if b == 0:
                    if r == 0:
                        pb_i[i] = 1.0
                    else:
                        pb_i[i] = 0.0
                if ((b <= 10.0) and (abs(i) <= 10)):
                    pb_i[i] = exp(-b*b)*special.iv(abs(r), b*b)
                else:
                    # otherwise, specia.iv() is probably numerically unstable
                    # instead compute this the ordinary way
                    pb_i[i] = exp(-0.5*(r/b)**2) / sqrt(2*pi*b*b)
                if r == 0:
                    pb_0 = pb_i[i]
            # Now take into account the "B" coefficient (for normalization).
            pb_i *= B / pb_0
            # Now add the two terms together
            p_i = pa_i - pb_i
            #plt.plot(r_i, p_i)
            plt.step(r_i+0.5, p_i)
            plt.show()

        elif sys.argv[1] == '-dogg':
            if i+6 >= len(sys.argv):
                raise InputError('Error: Expected 6 numeric arguments following \"'+sys.argv[1]+'\"\n')
            A = float(sys.argv[i+1])
            B = float(sys.argv[i+2])
            a = float(sys.argv[i+3])
            b = float(sys.argv[i+4])
            m = float(sys.argv[i+5])
            n = float(sys.argv[i+6])
            width = max(a,b)
            r_i = np.arange(-ceil(4.0*width), ceil(4.0*width), 1.0)
            p_i = np.array([A*exp(-(abs(r)/a)**m)-B*exp(-(abs(r)/b)**n) for r in r_i])
            #plt.plot(r_i, p_i)
            plt.step(r_i+0.5, p_i)
            plt.show()

        else:
            raise InputError('Error: You must select a filter type\n'
                             '       (eg "-gauss","-ggauss","-dog","-dogg")\n')
    except (InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return


if __name__ == '__main__':
    main(sys.argv)

