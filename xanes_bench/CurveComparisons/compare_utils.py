# This is a simple tool to compare to plots
# Fanchen Meng, 2022
# FC TODO 1. paser for different codes
#         2. compatible with convergece

import numpy as np
from scipy import interpolate
#from scipy.optimize import minimize_scalar, minimize
#import pathlib
#import os
#from os import environ as env
#from os import path
from scipy.stats import pearsonr
from scipy.stats import spearmanr
#import re

from scipy.ndimage import gaussian_filter1d

from xanes_bench.CurveComparisons.compare_pasers import XSplot_rescale

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

def comparePlots( omega, p1, p2, window=False, inverse=False):

    plot1, interp1 = makeInterpolateWithWindow( p1, window )
    plot2, interp2 = makeInterpolateWithWindow( p2, window )
    pearson = PearsonCoeff( omega, plot1, plot2, interp1, interp2, inverse )
    spearman = SpearmanCoeff( omega, plot1, plot2, interp1, interp2, inverse )
    coss = CosSimilar( omega, plot1, plot2, interp1, interp2, True)

    return pearson, spearman, coss
    #print( "{:12f}  {:12f}  {:10f}  {:12.8f}  {:12.8f}  {:12.8f}  {:10f}".format( radius, omega, alpha, coss, pearson, spearman, relArea, ))

def minCos(plot1, plot2, step=0.01):
    # hard coded
    # need to find a more robust way
    start = 12
    stop = -12

    if start <= stop:
        print('WARNING: Start {} is larger than stop {}]'.format(start, stop))
        exit()
    # 
    pearson={}
    spearman={}
    coss = {}

    i = start
    while i > stop:
        pearson[i],spearman[i],coss[i] = comparePlots(i, plot1, plot2)
        i -=step

    # find omega according to coss
    m = 100
    for i,j in coss.items():
        if j < m:
            m = j
            m_ind = i

    ## check if the gradient makes sense
    gplot1 = np.vstack((plot1[:,0],np.gradient(plot1[:,1],plot1[1,0]-plot1[0,0]))).T
    gplot2 = np.vstack((plot2[:,0],np.gradient(plot2[:,1],plot2[1,0]-plot2[0,0]))).T
    
    x1 = peak_loc(gplot1)
    x2 = peak_loc(gplot2)
    
    if abs(x1 - m_ind - x2 ) < 1:
        pass #return m_ind
    else:
        cossplot = np.vstack((list(coss.keys()),list(coss.values()))).T
        for i in valley_loc(cossplot):
            if abs(x1 - i - x2 ) < 1:
                print('*',end='.')
                m_ind = i

##     if m_ind == start :
##         print('{} Start wrong!'.format(mpid))
##     elif m_ind == stop:
##         print('{} Stop wrong!'.format(mpid))
        
    return pearson, spearman, coss, m_ind 

def peak_loc(plot):
    ''' locate the peak positon of a spectrum
        
        Parameters
        ----------
        plot : 2d-array

        Returns
        -------
        position of the peak
    '''
    return plot[plot[:,1].argmax(),0]

def valley_loc(plot):
    ''' locate the valley positions of a spectrum

        Parameters
        ----------
        plot : 2d-array


    '''
    valley = []
    valley_tmp = []
    length = int(np.ceil(len(plot)/60))
    for i in range(1,len(plot)-1):
        if plot[i,1] - plot[i-1,1] < 0 and plot[i+1,1] - plot[i,1] > 0:
            valley_tmp.append(i)

    for i in valley_tmp:
        if i < length:
            left = left = np.array(plot[0:i,1]).mean()
        else:
            left = np.array(plot[i-length:i,1]).mean()

        if i + length > len(plot):
            right = np.array(plot[i+1:len(plot),1]).mean()
        else:
            right = np.array(plot[i+1:i+length,1]).mean()


        if plot[i,1] < left and plot[i,1] < right:
                valley.append(plot[i,0])
    return valley

def trucate_spectra(plot):
    ''' truncate spectra from OCEAN or EXCITING

        Parameters
        ----------
        plot : one of the parsers defined in compare_parsers, mandatory
            spectra class

        Return
        ------
        index : int
    '''
    x = plot.nx
    y = plot.ny

    logic = y > 0.5

    seq = y == y[logic][-1]
    index = seq.argmax()
    return index


def compare_between_codes(file1, file2, itm=1):
    ''' Atomatic align the spectra and calculate the spearman coefficient
        Comparison was done for epsilon instead of cross-section

        Parameters
        ----------
        file1 : one of the parsers defined in compare_parsers, mandatory
            spectra1
        file2 : one of the parsers defined in compare_parsers, mandatory
            spectra2

        Return
        ------
        shift : real
            relative shift between the two spectra, sign is meaningful
        spearman : real
            spearman coefficient, which is calculated as log(1-spearman)
    '''
    
    if not isinstance(file1, XSplot_rescale):
        index = trucate_spectra(file1)
        plot1 = np.stack((file1.x[:index], file1.y[:index])).T
    else:
        plot1 = np.stack((file1.x, file1.y)).T

    if not isinstance(file2, XSplot_rescale):
        index = trucate_spectra(file2)
        plot2 = np.stack((file2.x[:index], file2.y[:index])).T
    else:
        plot2 = np.stack((file2.x, file2.y)).T

    pears, spear, coss, shift = minCos(plot1, plot2, step=0.01)
    
    spearman = comparePlots( shift, plot1, plot2)[itm] # 0 for pearson, 1 for spearman, 2 for coss
    if itm != 2:
        return shift, np.log10(1-spearman)
    else:
        return shift, np.log10(spearman)
