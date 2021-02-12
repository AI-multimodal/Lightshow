import json
import sys
import os
import pathlib
from os import environ as env
import re
from functools import reduce
import shutil
import subprocess
import math
import glob

#TODO: trim imports


def collectPW( rundir: str, savedir: str ):
    os.makedirs(savedir, exist_ok=True )
    os.chdir( rundir )
    # Set PW.x path

    
    # SCF run
    if not os.path.isfile("scf.in"):
        print( "Failed to find scf.in" )
        exit()
    if not os.path.isfile("scf.out"):
        print( "Failed to find scf.out" )
        exit()
    xmlFile = os.path.join( 'scf.save', 'data-file-schema.xml' )
    if not os.path.isfile( xmlFile ):
        print( "Failed to find {:s}".format( xmlFile ) )
        exit()

    targDir = pathlib.Path(savedir, 'scf')
    os.makedirs(targDir, exist_ok=True)
    shutil.copyfile( 'scf.in', os.path.join(targDir, 'scf.in') )
    shutil.copyfile( 'scf.out', os.path.join(targDir, 'scf.out') )
    shutil.copyfile( xmlFile, os.path.join(targDir, 'pwscf.xml' ) )
    

    #NSCF run
    if not os.path.isfile("nscf.in"):
        print( "Failed to find nscf.in" )
        exit()
    if not os.path.isfile("nscf.out"):
        print( "Failed to find nscf.out" )
        exit()
    xmlFile = os.path.join( 'nscf.save', 'data-file-schema.xml' )
    if not os.path.isfile( xmlFile ):
        print( "Failed to find {:s}".format( xmlFile ) )
        exit()

    targDir = pathlib.Path(savedir, 'nscf')
    os.makedirs(targDir, exist_ok=True)
    shutil.copyfile( 'nscf.in', os.path.join(targDir, 'nscf.in') )
    shutil.copyfile( 'nscf.out', os.path.join(targDir, 'nscf.out') )
    shutil.copyfile( xmlFile, os.path.join(targDir, 'pwscf.xml' ) )

    # DOS, but drop into NSCF
    for f in ["dos.in", "dos.out", "nscf.dos"]:
        if not os.path.isfile(f):
            print( "Failed to find {:s}".format(f) )
            exit()
        else:
            shutil.copyfile( f,  os.path.join(targDir, f ) )


    bandList = glob.glob("nscf_band*.in")

    for bandFile in bandList:
        outFile = bandFile.replace("in", "out" )
        if not os.path.isfile(bandFile):
            print( "Failed to find {:s}".format(bandFile) )
            exit()
        if not os.path.isfile(outFile):
            print( "Failed to find {:s}".format(outFile) )
            exit()
        xmlFile = os.path.join( bandFile.replace("in", "save"), 'data-file-schema.xml' )
        if not os.path.isfile( xmlFile ):
            print( "Failed to find {:s}".format( xmlFile ) )
            exit()

        targDir = pathlib.Path( savedir, bandFile.replace("nscf_","").replace(".","").replace("in","") )
        os.makedirs(targDir, exist_ok=True)
        shutil.copyfile( bandFile, os.path.join(targDir, bandFile) )
        shutil.copyfile( outFile, os.path.join(targDir, outFile) )
        shutil.copyfile( xmlFile, os.path.join(targDir, 'pwscf.xml' ) )



def main():

    mpid = "mp-" + input("Input the mp number: ")
    program = input("[O]cean, [X]spectra, or [E]xciting? ")
    method = 'g'
    print ("Method hardwired for ground state!")

    ### Need to find a home for this
    with open ('cluster.json', 'r') as fd:
        clusterJSON = json.load(fd)

    subdir = pathlib.Path(env['PWD'], "data", "mp_structures",mpid )
    if not os.path.isdir(subdir):
        print("Couldn't find ", subdir )
        exit()

    for p in program:
        if( p == 'O' or p == 'o' ):
            rundir = pathlib.Path(env['PWD'], "data", "mp_structures",mpid, 'OCEAN', 'groundState' )
            savedir = pathlib.Path(env['PWD'], "save", "mp_structures",mpid, 'OCEAN', 'groundState' )
            if not os.path.isdir(rundir):
                print("Couldn't find ", rundir )
                exit()
            collectPW( rundir, savedir )
        elif( p == 'X' or p == 'x' ):
            rundir = pathlib.Path(env['PWD'], "data", "mp_structures",mpid, 'XS', 'groundState' )
            savedir = pathlib.Path(env['PWD'], "save", "mp_structures",mpid, 'XS', 'groundState' )
            if not os.path.isdir(rundir):
                print("Couldn't find ", rundir )
                exit()
            collectPW( rundir, savedir )
        elif( p == 'E' or p == 'e' ):
            print("Exciting not implemented!")
            exit()
        else:
            print( "Didn't recognize program selected: ", program )


if __name__ == '__main__':
    main()
