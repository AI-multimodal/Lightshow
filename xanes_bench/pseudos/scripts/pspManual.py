import json
import os
import base64, bz2



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
  
    return True


def main():
    dbRootName = input( "What is the database root?")

    
    with open( dbRootName + '.json', 'r' ) as f:
        js = json.load( f )
    with open( dbRootName + '_pseudos.json', 'r' ) as f:
        pspjs = json.load( f )



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
    update = False
    if update:
        with open( dbRootName + '.json', 'w' ) as f:
            f.write( json.dumps( js ) )
        with open( dbRootName + '_pseudos.json', 'r' ) as f:
            f.write( json.dumps( pspjs  ))


if __name__ == '__main__':
  main()
