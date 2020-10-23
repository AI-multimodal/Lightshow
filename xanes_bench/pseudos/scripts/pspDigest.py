import json
import os
import base64, bz2




pspJsonIn = 'SSSP_precision.json'
pspDir = 'SSSP_precision_pseudos'
pspJsonOut = 'SSSP_precision_pseudos.json'

with open( pspJsonIn, 'r' ) as f:
    inJSON = json.load(f)

outJSON = dict()

for element in inJSON:
    print( element )
    
    pspFile = os.path.join(pspDir, inJSON[element]['filename'])
    with open( pspFile, 'r' ) as f:
        pspString = f.read()

#    outJSON[element] = base64.b64encode(zlib.compress(pspString.encode("utf-8"),9)).decode("utf-8")
    outJSON[inJSON[element]['filename']] = base64.b64encode(bz2.compress(pspString.encode("utf-8"),9)).decode("utf-8")


#print( outJSON)

with open( pspJsonOut, 'w' ) as f:
#    json.dump(outJSON, f )
    f.write( json.dumps( outJSON ) )
    



