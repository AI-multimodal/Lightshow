# This is a simple tool to compare to plots, currently hardwired for DOS

import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar
import pathlib
import os
from os import environ as env



def CosSimilar( omega, plot1, plot2, inverse=False ):
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=0.0 )
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=0.0 )

    norm1 = np.sqrt( np.dot( plot1[:,1], plot1[:,1] ) )
    norm2 = np.sqrt( np.dot( plot2[:,1], plot2[:,1] ) )

    plot1at2 = interp1( plot2[:,0]+omega )
    plot2at1 = interp2( plot1[:,0]-omega )

    cosSimilarityA = np.dot( plot1[:,1], plot2at1 ) 
    cosSimilarityB = np.dot( plot2[:,1], plot1at2 )

    cosSimilarity = np.sqrt( cosSimilarityA ) * np.sqrt( cosSimilarityB ) / ( norm1 * norm2 )

    if inverse:
      return 1 - cosSimilarity
    return cosSimilarity




def comparePlots( saveFiles, shiftAvg=True ):

    plot1 = np.loadtxt( saveFiles[0],  dtype=np.float64, usecols=(0,1))
    plot2 = np.loadtxt( saveFiles[1],  dtype=np.float64, usecols=(0,1))

    if shiftAvg:
        plot1Avg = np.mean( plot1[:,1] )
        plot1[:,1] = plot1[:,1] - plot1Avg
        plot2Avg = np.mean( plot2[:,1] )
        plot2[:,1] = plot2[:,1] - plot2Avg


    omega = 0.0
    cosSimilarity = CosSimilar( omega, plot1, plot2 )
    print( "Cosine similarity = {:f} ".format(cosSimilarity) )

    res = minimize_scalar( CosSimilar, args = ( plot1, plot2, True ), options={'xtol': 1e-8})
    if res.success:
        print( "Shift = {:f} eV".format(res.x) )
        print( "Cos S = {:f}".format(1-res.fun) )
    else:
        print( "Optmizing energy shift failed" )
        exit()



def main():

    mpid = "mp-" + input("Input the mp number: ")

    print("Hardwired to compare [O]cean and [X]spectra!")
    program = ['o', 'x' ]

    print( "Hardwired to compare the plot")
    pathTail = pathlib.Path('groundState', 'nscf', 'nscf.dos' )

    saveFiles = []
    for i in program:
        if( i == 'O' or i == 'o' ):
            saveFile = pathlib.Path(env['PWD'], "save", "mp_structures",mpid, 'OCEAN', pathTail )
            if not os.path.isfile(saveFile):
                print("Couldn't find ", saveFile )
                exit()
        elif( i == 'X' or i == 'x' ):
            saveFile = pathlib.Path(env['PWD'], "save", "mp_structures",mpid, 'XS', pathTail )
            if not os.path.isfile(saveFile):
                print("Couldn't find ", saveFile )
                exit()
        elif( i == 'E' or i == 'e' ):
            print("Exciting not implemented!")
            exit()
        else:
            print( "Didn't recognize program selected: ", program )

        saveFiles.append( saveFile )


#    print( str(saveFiles[0]), str(saveFiles[1]) )
    comparePlots( saveFiles )



if __name__ == '__main__':
    main()

