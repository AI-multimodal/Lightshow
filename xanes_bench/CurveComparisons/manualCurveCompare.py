# This is a simple tool to compare to plots

import sys
import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar, minimize
import pathlib
import os
from os import environ as env
from scipy.stats import pearsonr
from scipy.stats import spearmanr

from scipy.ndimage import gaussian_filter1d


def gaussAvgAndDiff( gWidth, omega, plot1, plot2 ):
    firstp1 = plot1[1,1]
    lastp1 = plot1[-1,1]
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp1,lastp1) )
    firstp2 = plot2[1,1]
    lastp2 = plot2[-1,1]
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp2,lastp2) )

    plot1at2 = interp1( plot2[:,0]+omega )
    plot2at1 = interp2( plot1[:,0]-omega )

    sigma = gWidth/(plot1[1,0]-plot1[0,0])
    convPlot1 = gaussian_filter1d(plot1[:,1], sigma, axis=0, order=0, output=None, mode='nearest', truncate=4.0)
    convPlot2 = gaussian_filter1d(plot2at1, sigma, order=0, output=None, mode='nearest', truncate=4.0)

    convPlot = (convPlot1+convPlot2)/2

    diff1 = plot1[:,1] - convPlot
    diff2 = plot2at1 - convPlot

    
    norm1 = np.sqrt( np.dot( diff1, diff1 ) )
    norm2 = np.sqrt( np.dot( diff2, diff2 ) )

    cosSimilarityA = np.dot( diff1, diff2 )

    cosSimilarity = cosSimilarityA /  ( norm1 * norm2 )
    print( cosSimilarity )



#    # TODO comment out
#    for i in range(len(convPlot1)):
#        print( plot1[i,0], plot1[i,1], convPlot1[i], plot2at1[i], convPlot2[i], convPlot[i], diff1[i], diff2[i])


    # TODO
    # 1. calculate cosine sim and pearson of the (plot1[:,0]-convPlot[:]) and (plot2at1-convPlot)
    # 2. Flip to 1 at 2 and repeat the sam3
    # 3. mean and report
    

def windowAverage( plot1, windowWidth, windowDelta ):
    firstp1 = plot1[1,1]
    lastp1 = plot1[-1,1]
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp1,lastp1) )

    windowMax = 0
    nWindow = round(windowWidth/windowDelta)
    for i in plot1[:,0]:
#        x = np.arange(i, nWindow, windowDelta)
#        x = np.arange(i, i+windowWidth, windowDelta )
        x = np.linspace( i, i+windowWidth, nWindow )
        y = interp1( x )
        total = 0
        for j in y:
            total += j
        total /= nWindow
        if total > windowMax:
            windowMax = total

#    print( windowMax )
    return windowMax

# Function to calculate the pearson 
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
def PearsonCoeffOld( omega, plot1, plot2, inverse=False ):
    firstp1 = plot1[1,1]
    lastp1 = plot1[-1,1]
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp1,lastp1) )
    firstp2 = plot2[1,1]
    lastp2 = plot2[-1,1]
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp2,lastp2) )

    norm1 = np.sqrt( np.dot( plot1[:,1], plot1[:,1] ) )
    norm2 = np.sqrt( np.dot( plot2[:,1], plot2[:,1] ) )

    plot1at2 = interp1( plot2[:,0]+omega )
    plot2at1 = interp2( plot1[:,0]-omega )

    p1 = pearsonr( plot2at1, plot1[:,1] )[0]
    p2 = pearsonr( plot1at2, plot2[:,1] )[0]

#    p = np.sqrt( p1 * p2 )
    p = 0.5 * ( p1 + p2 )

    if inverse:
      return 1 - p
    return p


# Function to calculate the pearson 
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
def PearsonCoeff( omega, plot1, plot2, interp1, interp2, inverse=False ):
#    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )
#    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )

    for start in range( len( plot1[:,0]) ):
        if plot1[start,0]-omega > plot2[0,0]:
            break
    for stop in reversed( range( len( plot1[:,0] ) ) ):
        if plot1[stop,0]-omega < plot2[-1,0]:
            break

    clipPlot1 = plot1[start:stop,1]
    plot2at1 = interp2( plot1[start:stop,0]-omega )

    for start in range( len( plot2[:,0]) ):
        if plot2[start,0]+omega > plot1[0,0]:
            break
    for stop in reversed( range( len( plot2[:,0] ) ) ):
        if plot2[stop,0]+omega < plot1[-1,0]:
            break

    clipPlot2 = plot2[start:stop,1]
    plot1at2 = interp1( plot2[start:stop,0]+omega )

    p1 = pearsonr( plot2at1, clipPlot1 )[0]
    p2 = pearsonr( plot1at2, clipPlot2 )[0]

    p = 0.5 * ( p1 + p2 )

    if inverse:
      return 1 - p
    return p


def SpearmanCoeff( omega, plot1, plot2, interp1, interp2, inverse=False ):
#    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )
#    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )
        
    for start in range( len( plot1[:,0]) ):
        if plot1[start,0]-omega > plot2[0,0]:
            break
    for stop in reversed( range( len( plot1[:,0] ) ) ):
        if plot1[stop,0]-omega < plot2[-1,0]:
            break
    
    clipPlot1 = plot1[start:stop,1]
    plot2at1 = interp2( plot1[start:stop,0]-omega )
    
    for start in range( len( plot2[:,0]) ):
        if plot2[start,0]+omega > plot1[0,0]:
            break
    for stop in reversed( range( len( plot2[:,0] ) ) ):
        if plot2[stop,0]+omega < plot1[-1,0]:
            break
    
    clipPlot2 = plot2[start:stop,1]
    plot1at2 = interp1( plot2[start:stop,0]+omega )

    p1 = spearmanr( plot2at1, clipPlot1 )[0]
    p2 = spearmanr( plot1at2, clipPlot2 )[0]
            
    p = 0.5 * ( p1 + p2 ) 
        
    if inverse:
      return 1 - p
    return p



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
def CosSimilarOld( omega, plot1, plot2, inverse=False ):
    firstp1 = 0.0 #plot1[1,1]
    lastp1 = plot1[-1,1]
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp1,lastp1) )
    firstp2 = 0.0 #plot2[1,1]
    lastp2 = plot2[-1,1]
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp2,lastp2) )

    norm1 = np.sqrt( np.dot( plot1[:,1], plot1[:,1] ) )
    norm2 = np.sqrt( np.dot( plot2[:,1], plot2[:,1] ) )

    plot1at2 = interp1( plot2[:,0]+omega )
    norm3 = np.sqrt( np.dot( plot1at2, plot1at2 ) )
    plot2at1 = interp2( plot1[:,0]-omega )
    norm4 = np.sqrt( np.dot( plot2at1, plot2at1 ) )

    cosSimilarityA = np.dot( plot1[:,1], plot2at1 ) 
    cosSimilarityB = np.dot( plot2[:,1], plot1at2 )


    print ( cosSimilarityA, norm1, norm4 )
    print ( cosSimilarityB, norm2, norm3 )

#    cosSimilarity = np.sqrt( cosSimilarityA ) * np.sqrt( cosSimilarityB ) / ( norm1 * norm2 )
    cosSimilarity = cosSimilarityA / ( 2.0 * norm1 * norm4 ) + cosSimilarityB / ( 2.0 * norm2 * norm3 )
    if inverse:
      return 1 - cosSimilarity
    return cosSimilarity



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
def CosSimilar( omega, plot1, plot2, interp1, interp2, inverse=False ):
#    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )
#    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )

    width1 = plot1[-1,0]-plot1[0,0]
    width2 = plot2[-1,0]-plot2[0,0]
    width = min( width1, width2 )

    # require 50% overlap?
    width *= 0.5
    
    for start in range( len( plot1[:,0]) ):
        if plot1[start,0]-omega > plot2[0,0]:
            break
    for stop in reversed( range( len( plot1[:,0] ) ) ):
        if plot1[stop,0]-omega < plot2[-1,0]:
            break

    clipPlot1 = plot1[start:stop,1]
    plot2at1 = interp2( plot1[start:stop,0]-omega )

    # Not sure about this one
    if plot1[stop,0]-plot1[start,0] < 0.005*width:
        print( plot1[stop,0], plot1[stop,0], width, omega)
        return 0

    for start in range( len( plot2[:,0]) ):
        if plot2[start,0]+omega > plot1[0,0]:
            break
    for stop in reversed( range( len( plot2[:,0] ) ) ):
        if plot2[stop,0]+omega < plot1[-1,0]:
            break
    
    clipPlot2 = plot2[start:stop,1]
    plot1at2 = interp1( plot2[start:stop,0]+omega )

    norm1 = np.sqrt( np.dot( clipPlot1, clipPlot1 ) )
    norm2 = np.sqrt( np.dot( clipPlot2, clipPlot2 ) )

    norm3 = np.sqrt( np.dot( plot1at2, plot1at2 ) )
    norm4 = np.sqrt( np.dot( plot2at1, plot2at1 ) )

    cosSimilarityA = np.dot( clipPlot1, plot2at1 )
    cosSimilarityB = np.dot( clipPlot2, plot1at2 )


#    print ( cosSimilarityA/(norm1*norm4), cosSimilarityA, norm1, norm4 )
#    print ( cosSimilarityB/(norm2*norm3), cosSimilarityB, norm2, norm3 )

    cosSimilarity = cosSimilarityA / ( 2.0 * norm1 * norm4 ) + cosSimilarityB / ( 2.0 * norm2 * norm3 )

    if plot2[stop,0]-plot2[start,0] < width:
#        print( cosSimilarity )
        cosSimilarity *= 0.9**(width/(plot2[stop,0]-plot2[start,0]) )
#        print( cosSimilarity )

    if inverse:
      return 1 - cosSimilarity
    return cosSimilarity

# using scaling alpha and shift omega, determine the rmsd between two curves
#TODO: normalize by the total length and sampling rate
def RMSDcurvesOld( alpha, omega, plot1, plot2 ):
    firstp1 = plot1[1,1]
    lastp1 = plot1[-1,1]
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp1,lastp1) )
    firstp2 = plot2[1,1]
    lastp2 = plot2[-1,1]
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp2,lastp2) )

    plot1at2 = interp1( plot2[:,0]+omega )
    plot2at1 = interp2( plot1[:,0]-omega )

    rmsd1 = 0.0
    for i in range(len( plot1[:,1] )):
        rmsd1 += ( plot1[i,1]-alpha*plot2at1[i] ) ** 2
    rmsd2 = 0.0
    for i in range(len( plot2[:,1] )):
        rmsd2 += ( alpha*plot2[i,1]-plot1at2[i] ) ** 2

    return ( 0.5 * ( np.sqrt( rmsd1 ) + np.sqrt( rmsd2 ) ) )


# using scaling alpha and shift omega, determine the rmsd between two curves
#TODO: normalize by the total length and sampling rate
def RMSDcurves( alpha, omega, plot1, plot2, interp1, interp2, Plot=False ):
#    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )
#    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )

    for start in range( len( plot1[:,0]) ):
        if plot1[start,0]-omega > plot2[0,0]:
            break
    for stop in reversed( range( len( plot1[:,0] ) ) ):
        if plot1[stop,0]-omega < plot2[-1,0]:
            break

    clipPlotRange1 = plot1[start:stop,0]
    clipPlot1 = plot1[start:stop,1]
    plot2at1 = interp2( plot1[start:stop,0]-omega )

    for start in range( len( plot2[:,0]) ):
        if plot2[start,0]+omega > plot1[0,0]:
            break
    for stop in reversed( range( len( plot2[:,0] ) ) ):
        if plot2[stop,0]+omega < plot1[-1,0]:
            break

    clipPlot2 = plot2[start:stop,1]
    plot1at2 = interp1( plot2[start:stop,0]+omega )

    rmsd1 = 0.0
    for i in range(len( clipPlot1 )):
        rmsd1 += ( clipPlot1[i]-alpha*plot2at1[i] ) ** 2
    rmsd2 = 0.0
    for i in range(len( clipPlot2 ) ):
        rmsd2 += ( alpha*clipPlot2[i]-plot1at2[i] ) ** 2

    if Plot:
        for i in range(len(clipPlot1)):
            print( clipPlotRange1[i], clipPlot1[i], plot2at1[i]*alpha )

    return ( 0.5 * ( np.sqrt( rmsd1 ) + np.sqrt( rmsd2 ) ) )



# using scaling alpha and shift omega, determine the rmsd between two curves
#TODO: normalize by the total length and sampling rate
def RMSDcurvesWindow( alpha, omega, plot1, plot2, interp1, interp2, Plot=False ):
#    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )
#    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
#                                    bounds_error=True )

    height1 = 0.1*windowAverage( plot1, 10, 0.1 )
    height2 = 0.1*windowAverage( plot2, 10, 0.1 )

    if Plot:
        print( height1, height2 )
#    return

    for windowStart in range( len( plot1[:,0]) ):
        if plot1[windowStart,1] > height1:
            break
    if plot1[-1,1] > height1:
        windowStop = len(plot1[:,0])-1
        if Plot:
            print( plot1[windowStart,0], plot1[windowStop,0]  )
    else:
        for windowStop in reversed( range( len( plot1[:,0] ) ) ):
            if plot1[windowStop,1] > height1:
                break
        if Plot:
            print( plot1[windowStart,0], plot1[windowStop,0], plot1[int( 0.8 * (windowStop-windowStart ) ) + windowStart,0] )
        windowStop = int( 0.8 * (windowStop-windowStart ) ) + windowStart

    for start in range( len( plot1[:,0]) ):
        if plot1[start,0]-omega > plot2[0,0]:
            break
    for stop in reversed( range( len( plot1[:,0] ) ) ):
        if plot1[stop,0]-omega < plot2[-1,0] and stop <= windowStop :
            break

    clipPlotRange1 = plot1[start:stop,0]
    clipPlot1 = plot1[start:stop,1]
    plot2at1 = interp2( plot1[start:stop,0]-omega )


    for windowStart in range( len( plot2[:,0]) ):
        if plot2[windowStart,1] > height2:
            break
    if plot2[-1,1] > height2:
        windowStop = len(plot2[:,0])-1
        if Plot:
            print( plot2[windowStart,0], plot2[windowStop,0]  )
    else:
        for windowStop in reversed( range( len( plot2[:,0] ) ) ):
            if plot2[windowStop,1] > height2:
                break
        if Plot:
            print( plot2[windowStart,0], plot2[windowStop,0], plot2[int( 0.8 * (windowStop-windowStart ) ) + windowStart,0] )
        windowStop = int( 0.8 * (windowStop-windowStart ) ) + windowStart

    for start in range( len( plot2[:,0]) ):
        if plot2[start,0]+omega > plot1[0,0]:
            break
    for stop in reversed( range( len( plot2[:,0] ) ) ):
        if plot2[stop,0]+omega < plot1[-1,0] and stop <= windowStop:
            break

    clipPlot2 = plot2[start:stop,1]
    plot1at2 = interp1( plot2[start:stop,0]+omega )

    rmsd1 = 0.0
    for i in range(len( clipPlot1 )):
        rmsd1 += ( clipPlot1[i]-alpha*plot2at1[i] ) ** 2
    rmsd2 = 0.0
    for i in range(len( clipPlot2 ) ):
        rmsd2 += ( alpha*clipPlot2[i]-plot1at2[i] ) ** 2

#    if Plot:
#        for i in range(len(clipPlot1)):
#            print( clipPlotRange1[i], clipPlot1[i], plot2at1[i]*alpha )

    return ( 0.5 * ( np.sqrt( rmsd1 ) + np.sqrt( rmsd2 ) ) )


def relAreaBetweenCurves( alpha, omega, plot1, plot2 ):
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=True )
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=True )

    for start in range( len( plot1[:,0]) ):
        if plot1[start,0]-omega > plot2[0,0]:
            break
    for stop in reversed( range( len( plot1[:,0] ) ) ):
        if plot1[stop,0]-omega < plot2[-1,0]:
            break

    clipPlot1 = plot1[start:stop,1]
    plot2at1 = interp2( plot1[start:stop,0]-omega )

    for start in range( len( plot2[:,0]) ):
        if plot2[start,0]+omega > plot1[0,0]:
            break
    for stop in reversed( range( len( plot2[:,0] ) ) ):
        if plot2[stop,0]+omega < plot1[-1,0]:
            break

    clipPlot2 = plot2[start:stop,1]
    plot1at2 = interp1( plot2[start:stop,0]+omega )

    area1 = 0
    norm1 = 0
    norm2 = 0

    for i in range(len( clipPlot1 )):
        area1 += abs( clipPlot1[i]-alpha*plot2at1[i] )
        norm1 += clipPlot1[i]
        norm2 += (alpha*plot2at1[i])
    
#    print( area1, norm1, norm2 )
    area1 = area1 / ( norm1 + norm2 )

    area2 = 0
    norm1 = 0
    norm2 = 0

    for i in range(len( clipPlot2 )):
        area2 += abs( alpha*clipPlot2[i]-plot1at2[i] )
        norm1 += alpha*clipPlot2[i]
        norm2 += plot1at2[i]

#    print( area2, norm1, norm2 )
    area2 = area2 / ( norm1 + norm2 )

    # No factor of 1/2 because we added the two norms above
    return( area1+area2)

# using scaling alpha and shift omega, determine the maximum distance between the two 
# curves at any given point
def RMSDsinglePointCurves( alpha, omega, plot1, plot2 ):
    firstp1 = plot1[1,1]
    lastp1 = plot1[-1,1]
    interp1 = interpolate.interp1d( plot1[:,0], plot1[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp1,lastp1) )
    firstp2 = plot2[1,1]
    lastp2 = plot2[-1,1]
    interp2 = interpolate.interp1d( plot2[:,0], plot2[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False, fill_value=(firstp2,lastp2) )

    plot1at2 = interp1( plot2[:,0]+omega )
    plot2at1 = interp2( plot1[:,0]-omega )

    rmsd1 = 0.0
    for i in range(len( plot1[:,1] )):
        r = ( plot2[i,1]-alpha*plot2at1[i] ) ** 2
        if r > rmsd1: 
            rmsd = r
    rmsd2 = 0.0
    for i in range(len( plot2[:,1] )):
        r = ( plot2[i,1]-alpha*plot1at2[i] ) ** 2
        if r > rmsd2:
            rmsd2 = r

    return ( 0.5 * ( np.sqrt( rmsd1 ) + np.sqrt( rmsd2 ) ) )



def makeInterpolateWithWindow( plot, doWindow ):

    heightScale = 0.1
    windowWidth = 10
    windowDelta = 0.1

    if not doWindow:
      return plot, interpolate.interp1d( plot[:,0], plot[:,1], assume_sorted=True, kind='cubic',
                                    bounds_error=False )

    height = heightScale*windowAverage( plot, windowWidth, windowDelta )


    for windowStart in range( len( plot[:,0]) ):
        if plot[windowStart,1] > height:
            break
    if plot[-1,1] > height:
        windowStop = len(plot[:,0])-1
    else:
        for windowStop in reversed( range( len( plot[:,0] ) ) ):
            if plot[windowStop,1] > height:
                break
        windowStop = int( 0.8 * (windowStop-windowStart ) ) + windowStart

    return plot[0:windowStop,:], interpolate.interp1d( plot[0:windowStop,0], plot[0:windowStop,1], 
                                                       assume_sorted=True, kind='cubic', bounds_error=False )

    

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
def comparePlots( p1, p2, window=False ):

    plot1, interp1 = makeInterpolateWithWindow( p1, window )
    plot2, interp2 = makeInterpolateWithWindow( p2, window )

    maxp1 = plot1[np.argmax( plot1[:,1] ), 0 ]
    maxp2 = plot2[np.argmax( plot2[:,1] ), 0 ]
    
    omega = maxp1-maxp2
    cosSimilarity = CosSimilar( omega, plot1, plot2, interp1, interp2 )
    print( "Cosine similarity = {:f} ".format(cosSimilarity) )

    bounds = [ [ plot1[0,0] -plot2[-1,0], plot1[-1,0]-plot2[0,0] ] ]
    res = minimize( CosSimilar, omega, args = ( plot1, plot2, interp1, interp2, True ), method='L-BFGS-B', bounds=bounds, options={'iprint':-1, 'eps': 1e-08 }, jac='3-point')
    if not res.success:
        res = minimize( CosSimilar, omega, args = ( plot1, plot2, interp1, interp2, True ), method='Powell', bounds=bounds )
    if res.success:
        print( "Shift = {:f} eV".format(res.x[0]) )
        print( "Cos S = {:f}".format(1-res.fun) )
    else:
        print( "Optmizing energy shift failed" )
        exit()

    omega = res.x[0]
    coss = 1-res.fun

    pearson = PearsonCoeff( omega, plot1, plot2, interp1, interp2 )
    spearman = SpearmanCoeff( omega, plot1, plot2, interp1, interp2 )

#    print( "Pearson = {:f}".format( PearsonCoeff( res.x[0], plot1, plot2 ) ) )
    print( "Pearson = {:f}".format( pearson ) )
#    print( "Spearman = {:f}".format(SpearmanCoeff( res.x[0], plot1, plot2 ) ) )
    print( "Spearman = {:f}".format( spearman ) )

#    windowAverage( plot1, 20, 0.2 )

#    print( RMSDcurves( 1.0, res.x[0], plot1, plot2 ) )

#    gaussAvgAndDiff( 2.0, res.x[0], plot1, plot2 )
    alpha = 1.0
    res2 = minimize( RMSDcurves, alpha, args = ( res.x[0], plot1, plot2, interp1, interp2 ), 
                     method='Nelder-Mead', options={'fatol': 1e-8})
    if res.success:
        print( "Alpha = {:e}".format(res2.x[0]) )
        print( "RMSD  = {:e}".format(res2.fun ) )

    alpha = res2.x[0]
    relArea = 100*relAreaBetweenCurves( alpha, omega, plot1, plot2 )
#    print( "Rel. area between = {:f} %".format( 100*relAreaBetweenCurves( alpha, omega, plot1, plot2 ) ) )
    print( "Rel. area between = {:f} %".format( relArea ) )

#    alpha = 1.0
    res2 = minimize( RMSDcurvesWindow, alpha, args = ( res.x[0], plot1, plot2, interp1, interp2 ), 
                     method='Nelder-Mead', options={'fatol': 1e-8})
    if res.success:
        print( "Alpha = {:e}".format(res2.x[0]) )
        print( "RMSD  = {:e}".format(res2.fun ) )

    alpha = res2.x[0]
    print( "Rel. area between = {:f} %".format( 100*relAreaBetweenCurves( alpha, omega, plot1, plot2 ) ) )

#    print( "{:f}  {:f}  {:f}  {:f}  {:f}  {:f}".format( omega, alpha, coss, pearson, spearman, relArea ))

    print( "{:f}  {:f}  {:f}  {:f}  {:f}  {:f}".format( omega, alpha, coss, pearson, spearman, relArea ))

    RMSDcurvesWindow( res2.x[0], res.x[0], plot1, plot2, interp1, interp2, True )

def main():


    saveFiles = []
    
    if len(sys.argv) == 3:
        saveFiles.append(sys.argv[1])
        saveFiles.append(sys.argv[2])
    else:
        saveFile = input("First file: ")
#    saveFile = 'ocean.1.1'
        saveFiles.append( saveFile )
        saveFile = input("Second file: ")
#    saveFile = 'xspectra.0.1'
        saveFiles.append( saveFile )

    plot1 = np.loadtxt( saveFiles[0],  dtype=np.float64, usecols=(0,1))
    plot2 = np.loadtxt( saveFiles[1],  dtype=np.float64, usecols=(0,1))
#    print( str(saveFiles[0]), str(saveFiles[1]) )
    comparePlots( plot1, plot2, True )



if __name__ == '__main__':
    main()

