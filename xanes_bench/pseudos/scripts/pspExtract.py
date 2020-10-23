import json
import os
import bz2, base64
import sys
import hashlib



element = sys.argv[1]
pspLibrary = 'SSSP_precision'
#pspJsonOut = 'SSSP_precision_pseudos.json'



pspJsonIn = pspLibrary + '.json'

with open( pspJsonIn, 'r' ) as f:
    inJSON = json.load(f)

if element not in inJSON:
    print( "Failed to find " + element )
    exit()

pspJsonIn = pspLibrary + '_pseudos.json'
print( inJSON[element]['filename'] )
fileName = inJSON[element]['filename']

with open( pspJsonIn, 'r' ) as f:
    pspJSON = json.load(f)

if fileName not in pspJSON:
    print( "Incomplete psp database" )
    exit()

pspString = bz2.decompress(base64.b64decode( pspJSON[fileName] ))

print( 'Expected hash:  ' + inJSON[element]['md5'] )
print( 'Resultant hash: ' + hashlib.md5( pspString ).hexdigest() )


with open( fileName, 'w' ) as f:
    f.write( pspString.decode("utf-8") )
    


