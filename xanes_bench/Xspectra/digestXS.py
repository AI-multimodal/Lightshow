import os
import pprint
import re
import copy
import json
import hashlib
import pathlib
from shutil import copy2


#TODO, this should take a unit cell and parse to find actual unique site list
def findSites( f, uc=None ):
    if uc is None:
        sites = []
        for s in os.listdir(f):
            if os.path.isdir( f / s ):
                sites.append( f / s)
        if len(sites) == 0:
            return None
        return sites
    else:
        return None

def parseES( sites ):
    if len(sites) == 0:
        print( "Empty sites!!")
        return None
    inputList = []
    for s in sites:
        newInput = {}
#        if not os.path.isfile( s /  "es.in" ):
#            return None
        try:
            fd = open( s /  "es.in", 'r' )
        except IOError:
            return None
        while True:
            line = fd.readline()

            if not line:
                break
            line = re.sub( "\!.*", "", line )
            m= re.search( "\&(\w+)", line )
            if m:
                key = m.group(1)
                newInput[key] = dict()
                while True:
                    line = fd.readline()
                    line = re.sub( "\!.*", "", line )
                    m= re.search( "^\s*\/", line )
                    if m:
                        break
                    m = re.search( "(\S+)\s+=\s+(\S+)", line )
                    if m:
                        newInput[key][m.group(1)] = m.group(2)
            else:
                m = re.search( "(\w+)\s+(\w+)", line )
                if m:
#                    print( m.group(1) )
                    if m.group(1) == 'K_POINTS':
                        line = fd.readline()
                        m = re.search( "(\d+\s+\d+\s+\d+)\s+(\d+\s\d+\s+\d+)", line )
                        if m:
                            newInput["kpoints"] = m.group(1)
                            newInput["kshift"] = m.group(2)
                    elif m.group(1) == 'CELL_PARAMETERS':
                        newInput["a1"] = fd.readline().rstrip()
                        newInput["a2"] = fd.readline().rstrip()
                        newInput["a3"] = fd.readline().rstrip()

        inputList.append( newInput )

    for i in range (len( inputList )-1):
        if not dictCompare( inputList[i], inputList[i+1] ):
            return None

    return inputList[0]



def parseDipoles( sites, photons ):
    inputList = []
    for s in sites:
        for d in photons:
            newInput = {}
            fd = open( s / d / 'xanes.in' , 'r' )
            while True:
                line = fd.readline()

                if not line:
                    break
                line = re.sub( "\!.*", "", line )
                m= re.search( "\&(\w+)", line )
                if m:
                    key = m.group(1)
                    newInput[key] = dict()
                    while True:
                        line = fd.readline()
                        line = re.sub( "\!.*", "", line )
                        m= re.search( "^\s*\/", line )
                        if m:
                            break
                        m = re.search( "(\S+)\s+=\s+(\S+)", line )
                        if m:
                            newInput[key][m.group(1)] = m.group(2)
                else:
                    m = re.search( "(\d+\s+\d+\s+\d+)\s+(\d+\s\d+\s+\d+)", line )
                    if m:
                        newInput["kpoints"] = m.group(1)
                        newInput["kshift"] = m.group(2)

        inputList.append( newInput )

    for i in range (len( inputList )-1):
        if not dictCompare( inputList[i], inputList[i+1], 'xepsilon' ):
            return None

    return inputList[0]


def dictCompare( d1, d2, exclude=None ):
    toplevel = list( set(d1.keys()).union(d2.keys()) )
    for i in toplevel:
        if exclude is not None:
            if re.search( exclude, i ):
                continue
        if i not in d1:
            print(i)
            return False
        elif i not in d2:
            print(i)
            return False
        elif type(d1[i]) is dict and type(d2[i]) is dict:
            if not dictCompare( d1[i], d2[i] ):
                return False
        elif d1[i] != d2[i]:
            print(i)
            return False

    return True

        


def unifiedInput( sites, photons ):

    FullInput = dict()
    FullInput['ES'] = parseES( sites )
    if FullInput['ES'] is None:
        return None
#    pp = pprint.PrettyPrinter(indent=4)
#    pp.pprint( FullInput['ES'] )

    FullInput['XS'] = parseDipoles( sites, photons )
    if FullInput['XS'] is None:
        return None

    jo = json.dumps(FullInput, indent = 4)  
    with open( 'unfiedInput', 'w' ) as fd:
         print(jo, file=fd)
    fd.close()
#    pp.pprint( FullInput['XS'] )

    return FullInput

def esExited( sites ):
    for s in sites:
        try:
            fd =  open ( s / 'es.out', 'r' )
        except:
            print( "Failed to open " + str(s) + '/es.out' )
            return False
        while True:
            line = fd.readline()
            if not line:
                return False
            m = re.search( 'JOB DONE', line )
            if m:
                break
        fd.close()
    return True


def esPseudoHashes( sites ):
    h = []
    for s in sites:
        t = []
        try:
            fd =  open ( s / 'es.out', 'r' )
        except:
            print( "Failed to open " + str(s) + '/es.out' )
            return False
        while True:
            line = fd.readline()
            if not line:
                break
            m = re.search( 'MD5 check sum:\s+(\S+)', line )
            if m:
                t.append( m.group(1) )
        fd.close()
        t = sorted( list(set( t )))
        if len( h ) == 0:
            for i in t:
                h.append(i)
        else:
            for i in t:
                if i not in h:
                    return None
    return h


def XSexited( sites, photons):
    for s in sites:
        for p in photons:
            try:
                fd = open( s /  p / 'xanes.out', 'r' )
            except:
                print( "Failed to open " + str(s) + '/' + p + '/xanes.out' )
                return False
            while True:
                line = fd.readline()
                if not line:
                    return False
                m = re.search( 'END JOB XSpectra', line )
                if m:
                    break
            fd.close()
            if not os.path.isfile( s / p / 'xanes.dat'):
                print( str(s) + '/' + p + '/xanes.dat not found!' )
                return False
    return True
                
# TODO: Need to make this return false if a file isn't there or fails to copy!
def copyFiles( sites, photons, f ):
    for s in sites:
        copy2( pathlib.Path(s / 'es.in'), f / 'es.in' )
        copy2( pathlib.Path(s / 'es.out'), f / 'es.out' )
        for p in photons:
            folder = f / p
            folder.mkdir(parents=True, exist_ok=True)
            copy2( pathlib.Path(s / p / 'xanes.in'), folder / 'xanes.in' )
            copy2( pathlib.Path(s / p / 'xanes.out'), folder / 'xanes.out' )
            copy2( pathlib.Path(s / p / 'xanes.dat'), folder / 'xanes.dat' )
        
        
    return True

def saveRun( localdir, targdir, photonDirs ):
#    targdir = '/Users/jtv1/Scratch/xanes_bench/Trash/'
#    photonDirs = ['dipole1', 'dipole2', 'dipole3' ]
    
    sites = findSites( localdir )
    if sites is None:
        return None
#    print(sites)

    ui = unifiedInput( sites, photonDirs )
    if ui is None:
        return None

    if not esExited( sites ):
        print( "es didn't all finish correctly")
        return None

    if not XSexited( sites, photonDirs ):
        print( "XS didn't finish all" )
        return None

    pseudoHashes = esPseudoHashes( sites )
    if pseudoHashes is None:
        print( "Failed to find pseudo hashes in es.out or they didn't match")
    ui['XS']['psp'] = pseudoHashes

#    print(json.dumps( ui, indent=2, sort_keys=True))
#    print( hashlib.sha1( json.dumps( ui, indent=2, sort_keys=True).encode() ).hexdigest() )
    outHash = hashlib.sha1( json.dumps( ui, indent=2, sort_keys=True).encode() ).hexdigest() 

#    folder =  pathlib.Path(targdir) / hashlib.sha1( json.dumps( ui, indent=2, sort_keys=True).encode() ).hexdigest() 
    folder =  pathlib.Path(targdir) / outHash
    folder.mkdir(parents=True, exist_ok=True)

    of = pathlib.Path( folder / "unifiedInputs.json" )
#    print( of )
    with open( of, 'w' ) as fd:
        fd.write( json.dumps( ui, indent=2, sort_keys=True) )
        fd.close()

    if not copyFiles( sites, photonDirs, folder ):
        print( "Copying failed, will skip!" )
        shutil.rmtree( folder )


    return outHash


def recursiveLoad( f,  targdir, photonDirs, i, m=2 ):
    if i > m:
        return
    else:
        i += 1
        oh = saveRun( f, targdir, photonDirs )
        if oh is not None:
            print( "Captured: " + oh )
        else:
            for s in os.listdir(f):
                if os.path.isdir( f / s ):
                    recursiveLoad( f / s, targdir, photonDirs, i, m )

def main():
    targdir = '/Users/jtv1/Scratch/xanes_bench/Trash/mp-390/XS/NON/'
    photonDirs = ['dipole1', 'dipole2', 'dipole3' ]
    f = pathlib.Path(os.environ["PWD"])

    i = 0
    recursiveLoad( f, targdir, photonDirs, i, 2 )
#    oh = saveRun( f, targdir, photonDirs )
#    if oh is not None:
#        print( "Captured: " + oh )
#    else:
#        for s in os.listdir(f):
#            if os.path.isdir( f / s ):
#                oh = saveRun( f / s, targdir, photonDirs )
#                if oh is not None:
#                    print( "Captured: " + oh )

if __name__ == '__main__':
    main()
