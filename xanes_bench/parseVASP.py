import numpy as np
import os

def readFermi(filename):
    with open(filename, 'r') as f:
        for line in f.readlines():
            if 'Fermi energy:' in line:
                print(line)
                Ef = line.split()[2]
    return float(Ef)

def readIBZ(filename):
    kvec_dict=dict()
    kvec_list=list()
    weight=dict()
    with open(filename, 'r') as f:
        f.readline()
        nk = int(f.readline())
        f.readline()
        for i in range(nk):
            line = f.readline()
            kx = float(line.split()[0])
            ky = float(line.split()[1])
            kz = float(line.split()[2])
            if (kx + kx) <1e-10:
                kx=0.0
            if (ky + ky) <1e-10:
                ky=0.0
            if (kz + kz) <1e-10:
                kz=0.0
            kvec = '{:16.8f} {:16.8f} {:16.8f}'.format(kx, ky, kz)
            kvec_list.append(kvec)
            kvec_dict[kvec] = np.array([kx, ky, kz])
            weight[kvec] = int(line.split()[3])


    return kvec_dict, kvec_list, weight

def readEIGEN(filename, kvec_list):

    eigenvals = {}
    with open(filename, 'r') as f:
        for _ in range(5):
            f.readline()
        nelectron, nk, nband= f.readline().split()

        nelectron = int(nelectron)
        nk = int(nk)
        nband = int(nband)
        for i in range(nk):
            f.readline()
            line = f.readline()
            kx = float(line.split()[0])
            ky = float(line.split()[1])
            kz = float(line.split()[2])
            if (kx + kx) <1e-10:
                kx=0.0
            if (ky + ky) <1e-10:
                ky=0.0
            if (kz + kz) <1e-10:
                kz=0.0
            # kvec = '{:16.8f} {:16.8f} {:16.8f}'.format(kx, ky, kz)
            kvec = kvec_list[i]
            # print()
            assert((float(kvec.split()[0]) - kx) < 1e-6)
            assert((float(kvec.split()[1]) - ky) < 1e-6)
            assert((float(kvec.split()[2]) - kz) < 1e-6)
            etmp = []
            for j in range(nband):
                etmp.append(float(f.readline().split()[1]))
            eigenvals[kvec] = np.array(etmp)
    occBands = int(nelectron/2)
    unoccBands = nband - occBands
    return eigenvals, occBands, unoccBands





