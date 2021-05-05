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

def factors(n):    
    return sorted(list(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))

def parsePWforParal( stdout: str, runText: str, maxMem: int, maxCPU: int ):
    poolEfficiency = 0.95
    intraPoolEfficiency = 0.95

    # Can't go wider than 1 process per band or everything breaks
    KSstates = int(re.search('number of Kohn-Sham states=\s+(\d+)', stdout ).group(1))
    if KSstates is None:
        print( "Failed to parse number of Kohn-Sham states in nscf_band.1.in run")
        exit()
    if KSstates < 1 or KSstates > 1000000:
        print( "Unrealistic number of Kohn-Sham states:", KSstate )
    NKpoints = int(re.search('number of k points=\s+(\d+)', stdout ).group(1))
    if NKpoints is None:
        print( "Failed to parse number of k points in nscf_band.1.in run")
        exit()
    if NKpoints < 1 or NKpoints > 1000000:
        print( "Unrealistic number of k points:", NKpoints )
    print( "{:s} has {:d} k-points and {:d} bands.".format( runText, NKpoints, KSstates ))
    #TODO memory checking to constrain the paralellization for spare mem
    r = re.search('Estimated max dynamical RAM per process >\s+(\d+\.\d+)\s+(\w+)', stdout)
    MemPerProcess = float(r.group(1))
    memUnits = r.group(2)
    if( memUnits == 'MB' ):
        MemPerProcess /= 1024
        memUnits = 'GB'

    r = re.search('wfc\\ \\(w\\.\\ buffer\\):\s+(\d+\\.\d+)\s+(\w+)',stdout)
    WfcMemPerProcess = float(r.group(1))
    memUnits = r.group(2)
    if( memUnits == 'MB' ):
        WfcMemPerProcess /= 1024
        memUnits = 'GB'

    kpointFactors = factors(NKpoints)
    npool = 1
    ncpus = maxCPU
    bestScore = 0.0

#   First see if the number of kpoints is exactly divisible by the MPI processes
    if maxCPU in kpointFactors:
        # Then check memory
        memGuess = ( MemPerProcess - WfcMemPerProcess ) \
                 + WfcMemPerProcess * ( poolEfficiency/maxCPU + (1-poolEfficiency) )
        if memGuess < maxMem:
            npool = maxCPU
            bestScore = 1.0
            print( "Mem guess: {:f} GB".format( memGuess ))

    if npool != maxCPU:
# loop over total cpus to use from max to half
# loop over cpu factors and calculate score
#   Count the time it takes the slowest pool (most kpoints)
#   small penalty for G parallel as opposed to k parallel
#   % penalty for leaving out CPUs entirely
#   Score is [ ceil(nkpoints/npools)*npools/nkpt]^-1 * ( 1-0.02*ncpus/npools) * ncpusUsed/ncpustot
# keep best, and score, exit loop over ncpus when used/tot drops below best score:w
        for tcpus in range(maxCPU,int(maxCPU/2),-1):
            if bestScore > (tcpus/maxCPU):
                break
            for tPool in factors(tcpus):
                # Will break if fewer than 1 band per processor in a pool
                if tcpus/tPool > KSstates:
                    continue
                if tPool > NKpoints:
                    continue
                tScore = NKpoints/(tPool*math.ceil(NKpoints/tPool)) \
                       * ( 1.0-0.05*tcpus/tPool ) * tcpus/maxCPU
                if tScore > bestScore:
                    memGuess = ( MemPerProcess - WfcMemPerProcess ) \
                             + WfcMemPerProcess * ( poolEfficiency/tPool + (1-poolEfficiency) )
                    memGuess *= ( intraPoolEfficiency*tPool/tcpus + (1-intraPoolEfficiency) )

                    if memGuess < maxMem:
                      bestScore = tScore
                      npool = tPool
                      ncpus = tcpus
                      print( "Mem guess: {:f} GB".format( memGuess ))
    return ncpus, npool;


def runPW( rundir: str, clusterJSON: dict ):

    os.chdir( rundir )
    # Set PW.x path
    pw = pathlib.Path( clusterJSON["QE_BIN"], "pw.x" )

    # SCF run
    if not os.path.isfile("scf.in"):
        print( "Failed to find scf.in" )
        exit()

    pathlib.Path("pwscf.EXIT").touch()
    print( "Test scf run")
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],"1", pw,"-inp","scf.in"],
                        capture_output=True,text=True)

    # Set up parallelization for SCF
    ncpus, npool = parsePWforParal( run.stdout, "SCF", clusterJSON["memory"], clusterJSON["MPI"] )

    # Run SCF calculation    
    print( "Run SCF with {:d} MPI tasks and {:d} pools".format( ncpus, npool ))
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],str(ncpus),
                          pw,"-inp","scf.in","-nk",str(npool)],
                         capture_output=True,text=True)

    # Write out the output file, and if non-zero, the errors
    with open( 'scf.out', 'w' ) as fd:
        fd.write(run.stdout)
    if len( run.stderr ) > 0 :
        with open( 'scf.err', 'w' ) as fd:
            fd.write(run.stderr)

    # Quit if the SCF didn't run to the end
    if not re.search('JOB DONE', run.stdout ):
        print( 'SCF did not complete correctly!' )
        exit()
    print("SCF completed successfully")

    # Move SCF run directory
    if os.path.isdir('scf.save'):
        shutil.rmtree('scf.save')
    shutil.move('pwscf.save', 'scf.save' )

    # Copy everything that isn't a wfc file
    if not os.path.exists('pwscf.save'):
        os.makedirs('pwscf.save')
    for f in os.listdir('scf.save'):
        if os.path.isfile(os.path.join('scf.save', f) ) and not re.search( "wfc\d+.dat", f ):
            shutil.copyfile( os.path.join('scf.save', f), os.path.join('pwscf.save', f) )

    pathlib.Path("pwscf.EXIT").touch()
    print( "Test nscf run")
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],"1",
                          pw,"-inp","nscf.in"],capture_output=True,text=True)

    # Set up parallelization for NSCF
    ncpus, npool = parsePWforParal( run.stdout, "NSCF", clusterJSON["memory"], clusterJSON["MPI"] )

    # re-copy data-file-schema.xml
    shutil.copyfile( os.path.join('scf.save', 'data-file-schema.xml'),
                     os.path.join('pwscf.save', 'data-file-schema.xml') )


    print( "Run NSCF with {:d} MPI tasks and {:d} pools".format( ncpus, npool ))
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],str(ncpus),
                          pw,"-inp","nscf.in","-nk",str(npool)],
                         capture_output=True,text=True)

    # Write out the output file, and if non-zero, the errors
    with open( 'nscf.out', 'w' ) as fd:
        fd.write(run.stdout)
    if len( run.stderr ) > 0 :
        with open( 'nscf.err', 'w' ) as fd:
            fd.write(run.stderr)

    # Quit if the SCF didn't run to the end
    if not re.search('JOB DONE', run.stdout ):
        print( 'non-SCF did not complete correctly!' )
        exit()
    print("NSCF completed successfully")


    # Move NSCF run directory
    if os.path.isdir('nscf.save'):
        shutil.rmtree('nscf.save')
    shutil.move('pwscf.save', 'nscf.save' )

            
    # TODO loop over k-point paths, but for now, just run the first
    bandList = sorted(glob.glob("nscf_band*.in"))
    print( "Run NSCF for {:d} band segments".format(len(bandList)))
    bandIter = 0

    for bandFile in bandList:
        bandIter += 1

        bandOutFile = bandFile.replace( "in", "out" )

        # Copy everything that isn't a wfc file
        if not os.path.exists('pwscf.save'):
            os.makedirs('pwscf.save')

        for f in os.listdir('scf.save'):
            if os.path.isfile(os.path.join('scf.save', f) ) and not re.search( "wfc\d+.dat", f ):
                shutil.copyfile( os.path.join('scf.save', f), os.path.join('pwscf.save', f) )


        pathlib.Path("pwscf.EXIT").touch()
        print( "Test NSCF band", bandFile)
        run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],"1", 
                              pw,"-inp",bandFile],
                            capture_output=True,text=True)

        # Set up parallelization for NSCF
        runString = "Band {:d}".format( bandIter )
        ncpus, npool = parsePWforParal( run.stdout, runString, clusterJSON["memory"], clusterJSON["MPI"] )

        # The test run messes up the xml, might also mess with Hubbard files 
        for f in os.listdir('scf.save'):
            if os.path.isfile(os.path.join('scf.save', f) ) and not re.search( "wfc\d+.dat", f ):
                shutil.copyfile( os.path.join('scf.save', f), os.path.join('pwscf.save', f) )
        print( "Run NSCF band {:d} with {:d} MPI tasks and {:d} pools".format( bandIter, ncpus, npool ))
        run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],str(ncpus),
                              pw,"-inp",bandFile,"-nk",str(npool)],
                             capture_output=True,text=True)
        # Write out the output file, and if non-zero, the errors
        with open( bandOutFile, 'w' ) as fd:
            fd.write(run.stdout)
        if len( run.stderr ) > 0 :
            with open( bandFile.replace( "in", "err" ), 'w' ) as fd:
                fd.write(run.stderr)

        # Quit if the SCF didn't run to the end
        if not re.search('JOB DONE', run.stdout ):
            print( 'non-SCF Band {:d} did not complete correctly!'.format( bandIter) )
            exit()

        # Move NSCF run directory
        bandDir = bandFile.replace( "in", "save" )
        if os.path.isdir(bandDir):
            shutil.rmtree(bandDir)
        shutil.move('pwscf.save', bandDir )    


    # Run DOS calculation    
    dos = pathlib.Path( clusterJSON["QE_BIN"], "dos.x" )
    print( "Run DOS with {:d} MPI tasks".format( clusterJSON["MPI"] ))
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],str(clusterJSON["MPI"]),
                          dos,"-inp","dos.in"], capture_output=True,text=True)

    # Write out the output file, and if non-zero, the errors
    with open( 'dos.out', 'w' ) as fd:
        fd.write(run.stdout)
    if len( run.stderr ) > 0 :
        with open( 'dos.err', 'w' ) as fd:
            fd.write(run.stderr)


  
        

def main():

    mpid =input("Input the mp number(s): ")
    program = input("[O]cean, [X]spectra, and/or [E]xciting? ")
    method = 'g'
    print ("Method hardwired for ground state!")

    ### Need to find a home for this
    with open ('cluster.json', 'r') as fd:
        clusterJSON = json.load(fd)

    for mp in mpid.split():
        subdir = pathlib.Path(env['PWD'], "data", "mp_structures", "mp-" + mp )
        print( subdir )
        if not os.path.isdir(subdir):
            print("Couldn't find ", subdir )
            exit()


    for m in mpid.split():
        mp = "mp-" + m
        print( "Running {:s}".format( mp ) )
        # Maybe make this more robust to keep someone from running 'ooooooooooo', but for now ...
        #  ie, could use a hash/dictionary
        for p in program:
            if( p == 'O' or p == 'o' ):
                print( "Running OCEAN" )
                rundir = pathlib.Path(env['PWD'], "data", "mp_structures",mp, 'OCEAN', 'groundState' )
                if not os.path.isdir(rundir):
                    print("Couldn't find ", rundir )
                    exit()
                runPW( rundir, clusterJSON )
                print( "OCEAN complete" )
                print( "######################" )
            elif( p == 'X' or p == 'x' ):
                print( "Running Xspectra" )
                rundir = pathlib.Path(env['PWD'], "data", "mp_structures",mp, 'XS', 'groundState' )
                if not os.path.isdir(rundir):
                    print("Couldn't find ", rundir )
                    exit()
                runPW( rundir, clusterJSON )
                print( "Xspectra complete" )
                print( "######################" )
            elif( p == 'E' or p == 'e' ):
                print("Exciting not implemented!")
    #            exit()
                print( "EXCITING complete" )
                print( "######################" )
            else:
                print( "Didn't recognize program selected: ", p )

        print( "Done with {:s}".format( mp ) )

if __name__ == '__main__':
    main()
