from collections import OrderedDict
import numpy as np
import re
import os
import json
import hashlib
import pathlib
from ase.units import Bohr


def dictCompare( d1, d2, exclude=None ):
    out = {}
    toplevel = list( set(d1.keys()).union(d2.keys()) )
    for i in toplevel:
        if exclude is not None:
            if re.search( exclude, i ):
                continue
        if i not in d1 or i not in d2:
#            print(i)
            out[i] = True
        elif type(d1[i]) is dict and type(d2[i]) is dict:
            o = dictCompare( d1[i], d2[i], exclude )
            if o:
                out[i] = o
        elif d1[i] != d2[i]:
#            print(i)
            out[i] = True

    return out


def grabFromDict( d1, dmap ):
    result = []
    for i in dmap:
        if isinstance( dmap[i], dict ):
            for j in grabFromDict( d1[i], dmap[i] ):
                result.append( j )
#            result.append( grabFromDict( d1[i], dmap[i] ) )
        else:
            result.append( d1[i] )
    return result

def headersFromDict( dmap ):
    result = []
    for i in dmap:
        if isinstance( dmap[i], dict ):
            for j in headersFromDict( dmap[i] ):
                result.append( i + '@' + j )
        else:
            result.append( i )
    return result

def addKden( js ):
    a1 = np.fromstring( js['ES']['a1'], count=3, sep=" ",dtype=np.float64 ) 
    a2 = np.fromstring( js['ES']['a2'], count=3, sep=" ",dtype=np.float64 ) 
    a3 = np.fromstring( js['ES']['a3'], count=3, sep=" ",dtype=np.float64 ) 

    vol = np.dot( np.cross( a1, a2 ), a3 )
    recip = [ np.linalg.norm( np.cross( a3, a2 ) / vol ),
              np.linalg.norm( np.cross( a3, a1 ) / vol ),
              np.linalg.norm( np.cross( a1, a2 ) / vol ) ]
    
    kpoints = np.fromstring( js['ES']['kpoints'], count=3, sep=" ",dtype=int )
    klen = kpoints[0]/recip[0]
    for i in [1,2]:
        if kpoints[i]/recip[i] < klen:
            klen = kpoints[i]/recip[i]

    js['ES']['klen'] = klen / Bohr

    kpoints = np.fromstring( js['XS']['kpoints'], count=3, sep=" ",dtype=int )
    klen = kpoints[0]/recip[0]
    for i in [1,2]:
        if kpoints[i]/recip[i] < klen:
            klen = kpoints[i]/recip[i]

    js['XS']['klen'] = klen / Bohr
    return

def loadUnifiedInput( d ):
    with open( d / 'unifiedInputs.json', 'r' ) as fd:
        js =  json.load( fd )
        addKden( js )
        return js
    return None

def main():
    f = pathlib.Path(os.environ["PWD"])
    runNames = []
    runs = []
    for d in os.listdir( f ):
        if os.path.isdir( f / d ):
            if re.search( '[0-9,a-f]{40}', d ):
#                print( d )
                j = loadUnifiedInput( f / d )
                if j is not None:
                    runs.append( j )
                    runNames.append(d)

    if( len( runs ) < 2 ):
        print( "1 or 0 runs were found!")
        exit()

    exclude = 'kpoints'
#    exclude = None
    out = {}
    for i in range(len(runs)):
        for j in range(i):
            out.update( dictCompare( runs[j], runs[i], exclude ) )
    out = OrderedDict( out )
#    print( out )
    

    outOptions = headersFromDict( out )
    colsKeep = [ 0 ]
    colSort = 0
    if len( outOptions ) == 0:
        print( "Found no differences!!!!" )
        exit()
    elif len( outOptions ) == 1:
        print( "Found 1 option" )
        print( headersFromDict( out ) )
        colsKeep = [0]
        colSort = 0
    else:
        print( "Found {:d} options".format( len( outOptions ) ) )
        print( headersFromDict( out ) )

    combinedHashAndCols = []
    for i in range(len(runs)):
        t = []
        s = grabFromDict( runs[i], out)
        for j in colsKeep:
            t.append( s[j] )
        t.append( runNames[i] )
        combinedHashAndCols.append( t )

#    for i in range(len(combinedHashAndCols)):
#      print( combinedHashAndCols[i])

    sorted_multi_list = sorted(combinedHashAndCols, key=lambda x: x[colSort])
    for i in range(len(sorted_multi_list)):
        print( i, sorted_multi_list[i] )


if __name__ == '__main__':
    main()
