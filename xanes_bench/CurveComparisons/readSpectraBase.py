# This is a simple tool to compare to plots

import numpy as np
import pathlib
from pathlib import Path
import os
from os import environ as env
from os import path
import re
import sys






#TODO: add parsing for EXCITING

# parsing the folder and spectra data from XSpectra
class XSplot():
    def __init__(self, folder, absorber=['0'], multiplicity=1):
        self.path = folder
        self.absorber = absorber
        self.multiplicity = multiplicity
    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            for polar in ['1', '2', '3']:
                # need to add a test to see if this path is a file
                path = self.path + "/" + iab + "/" + "dipole" + polar + "/xanes.dat"
                if Spectra is None:
                    Spectra = self.multiplicity * np.loadtxt( path, skiprows=4, usecols=(0,1)  )
                else:
                    Spectra[:,1] += self.multiplicity * np.loadtxt( path, skiprows=4, usecols=(1)  )
        return Spectra

# parsing the folder and spectra data from OCEAN
class OCEANplot():
    def __init__(self, folder, element, absorber=['0'], multiplicity=1):
        self.path = folder
        self.absorber = absorber
        self.element = element
        self.multiplicity = multiplicity

    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            print(iab)
            for polar in ['1', '2', '3']:
                p = int(polar)
                site = int(iab) + 1
                # need to test if this path is a file
                path = self.path + "/" + "absspct_{:s}.{:04d}_1s_{:02d}".format( self.element, site, p)
                if Spectra is None:
                    Spectra = self.multiplicity * np.loadtxt( path, skiprows=2, usecols=(0,2)  )
                else:
                    Spectra[:,1] += self.multiplicity * np.loadtxt( path, skiprows=2, usecols=(2)  )
        return Spectra


def str2list(a:str):
    '''
    take a string with the shape e.g. "['0','1','2']"
    convert it to a list -> ['0','1','2']
    works only with maximum of two digits
    '''

    out = []
    key = False
    for i in range(len(a)):
        if key:
            key = False
            continue
        if a[i].isnumeric():
            if i + 1 < len(a):
                if a[i+1].isnumeric():
                    out.append(a[i] + a[i+1])
                    key = True
                else:
                    out.append(a[i])
            else:
                out.append(a[i])
    return out


# below will not be used for now

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
        
