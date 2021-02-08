# J Vinson 2020
# C Vorwerk 2020
"""
Make the exciting input file and the photon files
"""

import os
import numpy as np
import json

def readEigval(filename):
    # first read only the global header to get number of states at each k-point
    # (nstsv) and the number of k-points (nkpt)
    nr=1
    block=1
    data=[]
    occdata=[]
    bdata=[]
    bdata2=[]
    kvec=[]
    for line in open(filename):
        if (nr ==1):
            nkpt=int(line.split()[0])
        elif (nr ==2):
            nstsv=int(line.split()[0])
        nr+=1
    # read the eigenvalues for each k-point 
    nr=1
    for line in open(filename):
        if (nr>3+(4+nstsv)*(block-1)+2) and (nr<3+(4+nstsv)*block-1):
            bdata.append(float(line.split()[1]))
            bdata2.append(float(line.split()[2]))
        elif (nr == 3+(4+nstsv)*(block-1)+1):
            kvec.append('{:16.8f} {:16.8f} {:16.8f}'.format(
                float(line.split()[3]),
                float(line.split()[2]), float(line.split()[1])))
        elif (nr == 3+(4+nstsv)*block):
            data.append(bdata)
            occdata.append(bdata2)
            bdata=[]
            block+=1
        nr+=1
    data.append(bdata)
    occdata.append(bdata2)

    data=np.asarray(data)
    occdata=np.asarray(occdata)
    
    # reshape the output to a directory
    out={}
    for i in range(data.shape[0]):
        out[kvec[i]]=data[i,:]
    
    # get number of occupied states
    nocc=0
    for i in range(occdata.shape[1]):
        if occdata[0,i] != 2.0:
            break
        nocc+=1
    unocc=nstsv-nocc
    return out, nocc, unocc


