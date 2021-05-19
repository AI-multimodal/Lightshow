import os
import pprint
import re
import copy
import json


#TODO, this should take a unit cell and parse to find actual unique site list
def findSites( uc=None ):
    if uc is None:
        sites = []
        for s in os.listdir("."):
            if os.path.isdir( s ):
                sites.append(s)
        return sites
    else:
        return None

def parseES( sites ):
    inputList = []
    for s in sites:
        newInput = {}
        fd = open( s + "/es.in", 'r' )
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



def parseDipoles( sites ):
    inputList = []
    for s in sites:
        for d in [1,2,3]:
            newInput = {}
            fd = open( s + "/dipole{:d}/xanes.in".format(d), 'r' )
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

        


def main():
    sites = findSites()

    FullInput = dict()
    FullInput['ES'] = parseES( sites )
    if FullInput['ES'] is None:
        exit()
    pp = pprint.PrettyPrinter(indent=4)
#    pp.pprint( FullInput['ES'] )

    FullInput['XS'] = parseDipoles( sites )
    if FullInput['XS'] is None:
        exit()

    jo = json.dumps(FullInput, indent = 4)  
    with open( 'unfiedInput', 'w' ) as fd:
         print(jo, file=fd)
    fd.close()
#    pp.pprint( FullInput['XS'] )


if __name__ == '__main__':
    main()
