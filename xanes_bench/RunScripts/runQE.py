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

def factors(n):    
    return sorted(list(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))


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

    # Can't go wider than 1 process per band or everything breaks
    KSstates = int(re.search('number of Kohn-Sham states=\s+(\d+)', run.stdout ).group(1))
    if KSstates is None:
        print( "Failed to parse number of Kohn-Sham states in nscf_band.1.in run")
        exit()
    if KSstates < 1 or KSstates > 1000000:
        print( "Unrealistic number of Kohn-Sham states:", KSstate )
    NKpoints = int(re.search('number of k points=\s+(\d+)', run.stdout ).group(1))
    if NKpoints is None:
        print( "Failed to parse number of k points in nscf_band.1.in run")
        exit()
    if NKpoints < 1 or NKpoints > 1000000:
        print( "Unrealistic number of k points:", NKpoints )
    print( "SCF has {:d} k-points and {:d} bands.".format( NKpoints, KSstates ))

    #TODO memory checking to constrain the paralellization for spare mem
    r = re.search('Estimated max dynamical RAM per process >\s+(\d+\.\d+)\s+(\w+)', run.stdout)
    MemPerProcess = float(r.group(1))
    memUnits = r.group(2)
    if( memUnits == 'MB' ):
        MemPerProcess /= 1024
        memUnits = 'GB'

    kpointFactors = factors(NKpoints)
    npool = 1
    ncpus = clusterJSON["MPI"]
    bestScore = 0.0
    if clusterJSON["MPI"] in kpointFactors:
        npool = clusterJSON["MPI"]
        bestScore = 1.0
        #TODO memcheck
    else:
# loop over total cpus to use from max to half
# loop over cpu factors and calculate score
#   Count the time it takes the slowest pool (most kpoints)
#   small penalty for G parallel as opposed to k parallel
#   % penalty for leaving out CPUs entirely
#   Score is [ ceil(nkpoints/npools)*npools/nkpt]^-1 * ( 1-0.02*ncpus/npools) * ncpusUsed/ncpustot
# keep best, and score, exit loop over ncpus when used/tot drops below best score:w
        for tcpus in range(clusterJSON["MPI"],int(clusterJSON["MPI"]/2),-1):
            if bestScore > (tcpus/clusterJSON["MPI"]):
                break
            for tPool in factors(tcpus):
                # Will break if fewer than 1 band per processor in a pool
                if tcpus/tPool > KSstates:
                    continue
                tScore = NKpoints/(tPool*math.ceil(NKpoints/tPool)) \
                       * ( 1.0-0.02*tcpus/tPool ) * tcpus/clusterJSON["MPI"]
                if tScore > bestScore:
                    bestScore = tScore
                    npool = tPool
                    ncpus = tcpus


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
    # TODO, with mem checking we'll need to re-run the test case
    # which will need to happen before we copy the density file over
    for f in os.listdir('scf.save'):
        if os.path.isfile(os.path.join('scf.save', f) ) and not re.search( "wfc\d+.dat", f ):
            shutil.copyfile( os.path.join('scf.save', f), os.path.join('pwscf.save', f) )

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

    # Copy everything that isn't a wfc file
    if not os.path.exists('pwscf.save'):
        os.makedirs('pwscf.save')
            
    # TODO loop over k-point paths, but for now, just run the first
    if not os.path.isfile("nscf_band.1.in"):
      print( "Failed to find nscf_band.1.in" )
      exit()
    for f in os.listdir('scf.save'):
        if os.path.isfile(os.path.join('scf.save', f) ) and not re.search( "wfc\d+.dat", f ):
            shutil.copyfile( os.path.join('scf.save', f), os.path.join('pwscf.save', f) )


    pathlib.Path("pwscf.EXIT").touch()
    print( "Test NSCF band 1")
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],"1", 
                          pw,"-inp","nscf_band.1.in"],
                        capture_output=True,text=True)
    # Can't go wider than 1 process per band or everything breaks
    KSstates = int(re.search('number of Kohn-Sham states=\s+(\d+)', run.stdout ).group(1))
    if KSstates is None:
        print( "Failed to parse number of Kohn-Sham states in nscf_band.1.in run")
        exit()
    if KSstates < 1 or KSstates > 1000000:
        print( "Unrealistic number of Kohn-Sham states:", KSstate )
    NKpoints = int(re.search('number of k points=\s+(\d+)', run.stdout ).group(1))
    if NKpoints is None:
        print( "Failed to parse number of k points in nscf_band.1.in run")
        exit()
    if NKpoints < 1 or NKpoints > 1000000:
        print( "Unrealistic number of k points:", NKpoints )
    print( "NSCF band 1 has {:d} k-points and {:d} bands.".format( NKpoints, KSstates ))

    #TODO memory checking to constrain the paralellization for spare mem
    r = re.search('Estimated max dynamical RAM per process >\s+(\d+\.\d+)\s+(\w+)', run.stdout)
    MemPerProcess = float(r.group(1))
    memUnits = r.group(2)
    if( memUnits == 'MB' ):
        MemPerProcess /= 1024
        memUnits = 'GB'

    kpointFactors = factors(NKpoints)
    npool = 1
    ncpus = clusterJSON["MPI"]
    bestScore = 0.0
    if clusterJSON["MPI"] in kpointFactors:
        npool = clusterJSON["MPI"]
        bestScore = 1.0
        #TODO memcheck
    else:
        for tcpus in range(clusterJSON["MPI"],int(clusterJSON["MPI"]/2),-1):
            if bestScore > (tcpus/clusterJSON["MPI"]):
                break
            for tPool in factors(tcpus):
                # Will break if fewer than 1 band per processor in a pool
                if tcpus/tPool > KSstates:
                    continue
                tScore = NKpoints/(tPool*math.ceil(NKpoints/tPool)) \
                       * ( 1.0-0.02*tcpus/tPool ) * tcpus/clusterJSON["MPI"]
                if tScore > bestScore:
                    bestScore = tScore
                    npool = tPool
                    ncpus = tcpus

    # The test run messes up the xml, might also mess with Hubbard files 
    for f in os.listdir('scf.save'):
        if os.path.isfile(os.path.join('scf.save', f) ) and not re.search( "wfc\d+.dat", f ):
            shutil.copyfile( os.path.join('scf.save', f), os.path.join('pwscf.save', f) )
    print( "Run NSCF band 1 with {:d} MPI tasks and {:d} pools".format( ncpus, npool ))
    run = subprocess.run([clusterJSON["para_prefix"], clusterJSON["task_flag"],str(ncpus),
                          pw,"-inp","nscf_band.1.in","-nk",str(npool)],
                         capture_output=True,text=True)
    # Write out the output file, and if non-zero, the errors
    with open( 'nscf_band.1.out', 'w' ) as fd:
        fd.write(run.stdout)
    if len( run.stderr ) > 0 :
        with open( 'nscf_band.1.err', 'w' ) as fd:
            fd.write(run.stderr)

    # Quit if the SCF didn't run to the end
    if not re.search('JOB DONE', run.stdout ):
        print( 'non-SCF Band 1 did not complete correctly!' )
        exit()

    # Move NSCF run directory
    if os.path.isdir('nscf_band.1.save'):
        shutil.rmtree('nscf_band.1.save')
    shutil.move('pwscf.save', 'nscf_band.1.save' )    
    

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


    if( program[0] == 'O' or program[0] == 'o' ):
        rundir = pathlib.Path(env['PWD'], "data", "mp_structures",mpid, 'OCEAN', 'groundState' )
        if not os.path.isdir(rundir):
            print("Couldn't find ", rundir )
            exit()
        runPW( rundir, clusterJSON )
    elif( program[0] == 'X' or program[0] == 'x' ):
        rundir = pathlib.Path(env['PWD'], "data", "mp_structures",mpid, 'XS', 'groundState' )
        if not os.path.isdir(rundir):
            print("Couldn't find ", rundir )
            exit()
        runPW( rundir, clusterJSON )
    elif( program[0] == 'E' or program[0] == 'e' ):
        print("Exciting not implemented!")
        exit()
    else:
        print( "Didn't recognize program selected: ", program )


if __name__ == '__main__':
    main()
