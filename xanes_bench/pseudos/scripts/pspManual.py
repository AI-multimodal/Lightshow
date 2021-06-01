import json
import os
import base64, bz2
import hashlib



def pspAdd( js, pspjs ):

    element = input( "Which element?" )
    if len(element) > 1:
        element = element[0].upper() + element[1].lower()
    else:
        element = element[0].upper()
#        element[1] = element[1].lower()

    
    if element in js:
        print( "{:s} exists in the database already!" )
        print( js[element] )
        overwrite = input( "Overwrite? [y]es/[n]o" ).lower()[0]
        if overwrite == 'n':
            return False

    Z = input( "What is Z val?")
    cutoff = input( "What is the cutoff?")
    rho = input( "What is the rho cutoff?")
    dual = float(rho)/float(cutoff)
    filename = input("What is the pseduo file name?")
    pspCollection = input("What is the reference for the pseudo?")
    saveWFC = input( "Save the wavefunction files too? (y/n) " )
    wfcShell = []
    wfcFile = []
    if saveWFC[0] == 'y' or saveWFC[0] == 'Y':
        while True:
            ws = input( "Which core hole? (1s, 2p) (any other entry to exit)" )
            if ws == '1s' or ws == '2p':
                wfcShell.append( ws )
                wfcFile.append( input("What is the wfc file name?") )

            else:
                break


    with open( filename, 'r' ) as f:
        pspString = f.read()

    filename = os.path.basename( filename )
    pspjs[filename] = base64.b64encode(bz2.compress(pspString.encode("utf-8"),9)).decode("utf-8")

    md5hash = hashlib.md5( pspString.encode("utf-8") ).hexdigest()




    if element not in js:
        js[element] = {}
    js[element]['Z_val'] = Z
    js[element]['cutoff'] = cutoff
    js[element]['dual'] = dual
    js[element]['filename'] = filename
    js[element]['md5'] = md5hash
    js[element]['pseudopotential'] = pspCollection
    js[element]['rho_cutoff'] = rho
  
    if len( wfcShell ) == 0:
        return True
    if 'wfc' not in js[element]:
        js[element]['wfc'] = {}
    for i in range(len( wfcShell)):
        with open( wfcFile[i], 'r' ) as f:
            wfcString = f.read()

        wfcFile[i] = os.path.basename( wfcFile[i] )
        pspjs[ filename + '_' + wfcShell[i] ] = base64.b64encode(bz2.compress(wfcString.encode("utf-8"),9)).decode("utf-8")
        wfchash = hashlib.md5( wfcString.encode("utf-8") ).hexdigest()
        js[element]['wfc'][wfcShell[i]] = { 'md5': wfchash, 'file': filename + '_' + wfcShell[i] }
    return True


def main():
    dbRootName = input( "What is the database root?")

    
    if os.path.isfile( dbRootName + '.json' ):
        with open( dbRootName + '.json', 'r' ) as f:
            js = json.load( f )
        with open( dbRootName + '_pseudos.json', 'r' ) as f:
            pspjs = json.load( f )
    else:
        js = {}
        pspjs = {}



    update = False
    c = 'c'
    while c != 'q':

        c = input( "[Q]uery, [A]dd, or [Q]uit?" ).lower()[0]

        if c == 'a':
            update = (update or pspAdd( js, pspjs ) )
        elif c == 'q':
            print( 'query' )

        c = 'q'

    print( update)

    if update:
        with open( dbRootName + '.json', 'w' ) as f:
            f.write( json.dumps( js ) )
        with open( dbRootName + '_pseudos.json', 'w' ) as f:
            f.write( json.dumps( pspjs  ))


if __name__ == '__main__':
  main()
