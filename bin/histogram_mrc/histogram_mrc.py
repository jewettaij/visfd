#!/usr/bin/env python
import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import mrcfile      #(if module not found, then run "pip install mrcfile")


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


def main():
    try:
        rescale01 = False
        nbins = -1
        mask_name = ''

        argv = [arg for arg in sys.argv]

        i = 1

        while i < len(argv):

            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')

            if argv[i] == '-n':
                if i + 1 >= len(argv):
                    raise InputError('Error: the \"' + argv[i] + '\" argument should be followed by a positive integer.\n')
                nbins = int(argv[i + 1])
                del argv[i:i + 2]

            elif (argv[i] == '-rescale'):
                rescale01 = True
                del argv[i:i + 1]
            elif (argv[i] == '-mask') or (argv[i] == '-m'):
                mask_name = argv[i+1]
                del argv[i:i + 2]
            else:
                i += 1

        if len(argv) > 2:
            raise InputError('Error: Two many arguments or unrecongized argument:\n'
                             '       '+str(argv[1:])+'\n')
        elif len(argv) < 2:
            raise InputError('Error: You must supply the name of a file in .MRC (.REC) format.\n')
                 
        file_name = argv[1]

        try:
            hdata = []
            sys.stderr.write('Reading MRC file "'+file_name+'"\n')
            with mrcfile.open(file_name, 'r') as mrc:
                #sys.stdout.write('mrc.data.shape = '+str(mrc.data.shape)+'\n')
                hdata = mrc.data.flatten()
                if mask_name != '':
                    sys.stderr.write('Reading MRC file "'+mask_name+'" (mask)\n')
                    with mrcfile.open(mask_name, 'r') as mask:
                        if mask.data.shape != mrc.data.shape:
                            raise InputError('Error: The MRC files ("'+
                                             file_name+'" and "'+mask_name+'")\n'
                                             '       must have the same number of voxels in the x,y,z directions.\n'
                                             '       (Currently: '+str(mrc.data.shape)+' and '+str(mask.data.shape)+', respectively)\n')
                        mdata = map(mask.data.flatten())
                        # Each entry in mdata is assumed to be either 0 or 1
                        # (or rarely, a number between 0 and 1).
                        # We can restrict the entries to entries in the original
                        # hdata array whose corresponding mask value is non-zero
                        # by multiplying each entry in the two arrays together.
                        hdata *= mdata
                hmin = min(hdata)
                hmax = max(hdata)
                
                if hmin == hmax:
                    raise InputError('Error: The image has only one intensity value: '+str(hmin)+'\n')
                if rescale01:
                    # change to type float
                    hdata = np.array(hdata, dtype='float32')
                    hdata -= float(hmin)
                    hdata *= 1.0/(float(hmax)-float(hmin))

                if nbins < 0:
                    #If unspecified, guess a reasonable number of histogram bins
                    nbins = 32 # default (for 8 bit integers)
                    if ((mrc.data.dtype == np.int8) or   # 8 bit integers?
                        (mrc.data.dtype == np.uint8)):   # 8 bit integers?
                        nbins = (1 + int(hmax) - int(hmin))

                if rescale01:
                    hmin = 0.0
                    hmax = 1.0

                sys.stderr.write('nbins = ' + str(nbins) + '\n')
                sys.stderr.write('(hmin,hmax) = ('+str(hmin)+','+str(hmax)+')\n')
                delta_bin = (float(hmax) - float(hmin)) / float(nbins-1)
                sys.stderr.write('delta_bin = ' + str(delta_bin) + '\n')
                bins = [hmin + (i-0.5)*delta_bin for i in range(0, nbins+1)]
                assert(len(bins) >= 2)
                #sys.stderr.write('bins = ' + str(bins) + '\n')

            plt.hist(hdata,
                     bins,
                     normed=1,
                     facecolor='green',
                     alpha=0.75)
            plt.ylabel('Frequency')
            plt.xlabel('Intensity (note: Dark voxels in IMOD have low intensity.)')
            plt.grid(True)
            plt.show()

        except IOError:
            raise InputError('Error: Unable to open file \"' + file_name + '\" for reading.\n')
    except (InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return


if __name__ == '__main__':
    main()
