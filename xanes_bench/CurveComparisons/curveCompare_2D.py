# This is a simple tool to compare to plots

import numpy as np
from scipy import interpolate
from scipy.optimize import minimize_scalar, minimize
import pathlib
import os
from os import environ as env
from os import path
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import re
import sys

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


def relAreaBetweenCurves( alpha, omega, plot1, plot2, interp1, interp2 ):
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
def comparePlots( p1, p2, window=False, coss=[], pearson=[], spearman=[], relArea=[], inverse=False):

    plot1, interp1 = makeInterpolateWithWindow( p1, window )
    plot2, interp2 = makeInterpolateWithWindow( p2, window )

    maxp1 = plot1[np.argmax( plot1[:,1] ), 0 ]
    maxp2 = plot2[np.argmax( plot2[:,1] ), 0 ]
    
    omega = maxp1-maxp2

    bounds = [ [ plot1[0,0] -plot2[-1,0], plot1[-1,0]-plot2[0,0] ] ]
    res = minimize( CosSimilar, omega, args = ( plot1, plot2, interp1, interp2, True ), 
                    method='L-BFGS-B', bounds=bounds, options={'iprint':-1, 'eps': 1e-08 }, jac='3-point')
    if not res.success:
        res = minimize( CosSimilar, omega, args = ( plot1, plot2, interp1, interp2, True ), 
                        method='Powell', bounds=bounds )
    if not res.success:
        print( "Optmizing energy shift failed" )
        exit()

    omega = res.x[0]
    if inverse:
        coss.append( res.fun )
    else:
        coss.append( 1-res.fun )


    pearson.append( PearsonCoeff( omega, plot1, plot2, interp1, interp2, inverse ) )
    spearman.append( SpearmanCoeff( omega, plot1, plot2, interp1, interp2, inverse ) )

    alpha = 1.0
    res2 = minimize( RMSDcurves, alpha, args = ( res.x[0], plot1, plot2, interp1, interp2 ), 
                     method='Nelder-Mead', options={'fatol': 1e-8})
    if not res2.success:
        print( "Optimizing RMSD scaling failed")
        exit()

    alpha = res2.x[0]
    rmsd = res2.fun
    relArea.append( relAreaBetweenCurves( alpha, omega, plot1, plot2, interp1, interp2 ) )

def parseEXCITINGFile( string, polar ):
    plot = None
    for p in polar:
        f = "EPSILON-{:s}/EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC{:d}{:d}.OUT".format( string, p, p )
        if os.path.isfile( f ):
            if plot is None:
                plot = np.loadtxt( f, skiprows=18, usecols=(0,2) )
            else:
                plot[:,1] += np.loadtxt( f, skiprows=18, usecols=(2) )
        else:
            return None

    return plot


def parseOCEANFile( site, polar, kpoint1, kpoint2=None ):
    plot = None
    k1 = [ int(kpoint1[0]), int(kpoint1[1]), int(kpoint1[2]) ]
    if kpoint2 is not None:
        k2 = [ int(kpoint2[0]), int(kpoint2[1]), int(kpoint2[2]) ]
    for s in site:
        for p in polar:
            if kpoint2 is not None:
                f = "NGKPT-{:d}.{:d}.{:d}/CNBSE-{:d}.{:d}.{:d}/absspct_Ti.{:04d}_1s_{:02d}".format( 
                    k1[0], k1[1], k1[2], k2[0], k2[1], k2[2], s, p )
            else:
                f = "CNBSE-{:d}.{:d}.{:d}/absspct_Ti.{:04d}_1s_{:02d}".format( k1[0], k1[1], k1[2], s, p )
            if os.path.isfile( f ):
                if plot is None:
                    plot = np.loadtxt( f, skiprows=2, usecols=(0,2),  dtype=np.double )
                else:
                    plot[:,1] += np.loadtxt( f, skiprows=2, usecols=(2), dtype=np.double  )
            else:
                return None

    if kpoint2 is not None:
        f = "NGKPT-{:d}.{:d}.{:d}/CNBSE-{:d}.{:d}.{:d}/avg.txt".format(
                    k1[0], k1[1], k1[2], k2[0], k2[1], k2[2] )
    else:
        f = "CNBSE-{:d}.{:d}.{:d}/avg.txt".format( k1[0], k1[1], k1[2] )
    with open(f,"w") as f:
        for i in range(len(plot[:,1])):
            f.write( "{:e} {:e}\n".format(plot[i,0], plot[i,1]) )
    return plot

def parseXSpectraFile( site, polar, kpoint1, kpoint2 ):
    plot = None
    k1 = [ int(kpoint1[0]), int(kpoint1[1]), int(kpoint1[2]) ]
    k2 = [ int(kpoint2[0]), int(kpoint2[1]), int(kpoint2[2]) ]
    #TODO site support!!
    # Need to build out support for *correctly* averaging over symmetry unique sites
    # For now, 
    for s in site:
        for p in polar:
            f = "k-{:d}-{:d}-{:d}/Spectra-{:d}-{:d}-{:d}/{:d}/dipole{:d}/xanes.dat".format( 
                k1[0], k1[1], k1[2], k2[0], k2[1], k2[2], s, p )
            if os.path.isfile( f ):
                if plot is None:
                    plot = np.loadtxt( f, skiprows=4, usecols=(0,1) )
                else:
                    plot[:,1] += np.loadtxt( f, skiprows=4, usecols=(1) )
            else:
                return None

    return plot
        

def parseFileK( program, site, polarization, kpoint1, kpoint2 ):
    if program == 'O':
        return parseOCEANFile( site, polarization, kpoint1, kpoint2 ) 
    elif program == 'X':
        return parseXSpectraFile( site, polarization, kpoint1, kpoint2 )
    else:
        return None

def loadPlotsKpoints( program, site, polarization ):

    plots = []
    foundKpoints1 = []
    foundKpoints2 = []

    kpointArray = np.loadtxt( './k.txt', skiprows=1, usecols=(0,1,2,5) )

    for k1 in kpointArray:
        for k2 in kpointArray:
            p = parseFileK( program, site, polarization, k1[0:3], k2[0:3] )
            if p is not None:
                plots.append( p )
                foundKpoints1.append( tuple(k1) )
                if tuple(k2) not in foundKpoints2:
                    foundKpoints2.append( tuple(k2) )  


    return plots, foundKpoints1, foundKpoints2


# Given a string, find all the 
def loadPlotsS( program, site, polarization, string ):
    plots = []
    foundParam = []
    attemptParam = []

    if program == 'E':
        for d in os.listdir():
            t = re.search( "EPSILON-" + string + "(\d+\.?\.*)", d )
            if t:
                attemptParam.append( t.group(1) )

        for i in sorted( attemptParam ):
            p = parseEXCITINGFile( string+i, polarization )
            if p is not None:
                plots.append( p )
                # This hack is because the k-points have the length param as the 4th object in a list
                foundParam.append( [0,0,0,float(i)])

    return plots, foundParam

# Print the similarity metrics to files
def printdata(coss, pearson, spearman, relArea, kpoints1, kpoints2, plots, inverse=False, site = 0):
    if inverse:
        coss = np.log10(coss)
        pearson = np.log10(pearson)
        spearman = np.log10(spearman)
        relArea = np.log10(relArea)
    lookup = {"COS.dat":coss, "Pearson.dat":pearson, "Spearman.dat":spearman, "relArea.dat":relArea}
    for file_name in lookup:
        filename = file_name.split(".")[0]+str(site[0])+".dat"
        with open(filename,"w") as f:
            print("    " + filename + "    ", file=f, end=",")
            for i in range(len(kpoints2)):
                print("{:12f}".format( kpoints2[i][3] ), file=f, end=",")	
            for i in range(len(plots)):
                if i % len(kpoints2) == 0 :
                    print("\n {:12f}".format( kpoints1[i][3] ), file=f, end = ",")
                    print("{:12.8f}".format( lookup[file_name][i] ), file=f, end="," )
                else:
                    print("{:12.8f}".format( lookup[file_name][i] ), file=f, end="," )
            f.write("\n")

def main():
    if len(sys.argv) > 1:
        inverse = True
    else:
        inverse = False
    program = 'X'
    method = 'K'
    string = 'gk'
#    site = [1,2,3,4,5,6,7,8]
    site = [0]
    polarization = [1,2,3]

    if method == 'K':
        plots, kpoints1, kpoints2 = loadPlotsKpoints( program, site, polarization )
    elif method == 'S':
        plots, kpoints = loadPlotsS( program, site, polarization, string )

    if len(plots) < 2:
        print( "Didn't find enough plots" )
        exit()

    print( "#i Largest k-point grid in es.in:    ",  kpoints1[-1])
    print( "#i Largest k-point grid in xanes.in: ",  kpoints2[-1])

    # Radius is replaced with the param for non-kpoint runs
    # print( "#   Radius       Shift        Scale       CosSimilar  Pearson     Spearman   Rel Area" )
    # print( "#   (Bohr)       (eV)                                                          (%)" )
    #print(kpoints2[:])
#    print("    SPEARMAN     ", end=" ")
#    for i in range(len(kpoints2)):
#        print("{:12f}".format( kpoints2[i][3] ), end=" ")

#    for i in range(len(plots)):
#        if i % len(kpoints2) == 0 :
#            swtch = True
#        else:
#            swtch = False 
    coss = []
    pearson = []
    spearman = []
    relArea= []
    for i in range(len(plots)):
        comparePlots( plots[-1], plots[i], True, coss, pearson, spearman, relArea, inverse)
    printdata(np.array(coss), np.array(pearson), np.array(spearman), np.array(relArea), kpoints1, kpoints2, plots, inverse, site)
    exit()

#    if swtch:
#        print("\n {:12f}".format( radius ), end = " ")
#        print("{:12.8f}".format(spearman), end=" " )
#    else:
#        print("{:12.8f}".format(spearman), end=" " )
#



if __name__ == '__main__':
    main()

