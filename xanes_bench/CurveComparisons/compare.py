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
from curveCompareBase import comparePlots, loadPlotsKpoints, loadPlotsS, printdata
from readSpectraBase import buildOCEANPath, buildXSpectraPath, parseOCEANFile, parseXSpectraFile


def main():
    if len(sys.argv) < 5 :
        print( "Usage: Path-1 Path-2 Site Polarization")
        exit()

    paths = [ sys.argv[1], sys.argv[2] ]
    site = int( sys.argv[3] )
    polar = int( sys.argv[4] )

    plots = []
    for i in range(2):
        f = buildOCEANPath( paths[i], 'Ti', site, polar )
        if( f.is_file() ):
            print( f )
            plots.append( parseOCEANFile( f ) )
            continue
        f = buildXSpectraPath( paths[i], site-1, polar )
        if( f.is_file() ):
            print( f )
            plots.append( parseXSpectraFile( f ) )
            continue

        print("Failed!")


    coss = []
    pearson = []
    spearman = []
    relArea= []
    alpha=[]
    omega=[]

    comparePlots( plots[0], plots[1], True, coss, pearson, spearman, relArea, omega, alpha )

    print( omega[0], alpha[0] )
    print( coss[0] )
    print( pearson[0] )
    print( spearman[0] )
    print( relArea[0] )




if __name__ == '__main__':
    main()

