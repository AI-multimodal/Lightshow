# This is a simple tool to compare to plots, currently hardwired for DOS

import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar
import pathlib
import os
from os import environ as env


# Function to calculate the cosine similarity between two curves (vectors)
# omega: relative energy shift between the two curves
# plot1: curve one as a numpy array, plot1[:,0] are the energies, plot1[:,1] are the values
# plot2: same as plot1
# inverse: boolean to return (1 - cosine similarity) instead (useful for optimization)
#
# plot1 and plot2 don't need to be on the same energy grid, a cubic interpolation is used
#
# LIMITATIONS: 
# 1. Assumes that the plots are on uniform energy grids, otherwise the cossimilarity should include a dx(?)
#
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



# Takes a list of two files, assumes both files are formated with two columns 
# compares the two datasets using cosineSimilarity
#
# saveFiles: list of files to read
# shiftAvg: Option to shift the y-values of each data set such that the average value is 0
#
#TODO: 
#  1. Column flexibility, e.g., 1:3 [(0,2)]
#  2. Different comparison methods
#  3. Default shift for the search
#
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

