# This is a simple tool to compare to plots

import numpy as np
import pathlib
from pathlib import Path
import os
from os import environ as env
from os import path
import re
import sys




#def buildEXCITINGPath( root, 

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

# Returns a path object pointing to an OCEAN spectra file for element, site, and polarization 
#  living in the root directory
def buildOCEANPath( root, element, site, polar ):
    f = "absspct_{:s}.{:04d}_1s_{:02d}".format( element, site, polar )
    return Path( root ) / f


# Parses an ocean spectrum and returns
# f is a Path object to the file
# if plot is present sum
# multiplicity is for the future when we support reduction in site/polarization symmetry
def parseOCEANFile( f, plot=None, multiplicity=1 ):
    if not f.is_file():
        return False

    if plot is None:
        plot = np.loadtxt( f, skiprows=2, usecols=(0,2),  dtype=np.double )
        if multiplicity != 1:
            plot[:,1] *= multiplicity
    else:
        plot[:,1] += multiplicity * np.loadtxt( f, skiprows=2, usecols=(2), dtype=np.double  )

    return plot


#
def buildXSpectraPath( root, site, polar ):
    s = "{:d}".format( site )
    d = "dipole{:d}".format(polar)
    return Path( root ) / s / d / "xanes.dat"


def parseXSpectraFile( f, plot=None, multiplicty=1 ):
    if not f.is_file(): 
        return False

    if plot is None:
        plot = np.loadtxt( f, skiprows=4, usecols=(0,1) )
    else:
        plot[:,1] += np.loadtxt( f, skiprows=4, usecols=(1) )

    return plot
        
