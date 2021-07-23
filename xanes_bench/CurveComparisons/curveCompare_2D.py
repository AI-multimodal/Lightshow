import os
import sys
import numpy as np
import pandas as pd
from readSpectraBase import XSplot, OCEANplot, str2list
from curveCompareBase import comparePlots

def find_XS_folder(folder):
    '''
    used to find the folders for ki and kf for the 2D comparison
    return: gs, es, absorber
    '''
    def find_gs_folder():
        gslst = []
        for i in os.listdir(folder):
            if i.startswith('k-'):
                gslst.append(i)
        return gslst

    def find_es_folder(gslst):
        eslst = []
        for i in gslst:
            path = folder + "/" + i
            for j in os.listdir(path):
                if j.startswith('Spectra-'):
                    eslst.append(j)
        eslst = list(set(eslst))
        return eslst

    def find_absorber(gslst, eslst):
        absorber = []
        for i in eslst:
            path = folder + "/" + gslst[0] + "/" + i
        if os.path.isdir(path):
            for j in os.listdir(path):
                if j.isnumeric():
                    absorber.append(j)
        else:
            eslst.remove(i)
        absorber = list(set(absorber))
        return eslst, absorber

    def screen_folder(gslst, eslst, absorber):
        for i in gslst:
            criteria1 = []

            for j in eslst:
                #criteria2 = []
                for k in absorber:
                    for polar in [1,2,3]:
                        path = folder + "/" + i + "/" + j + "/" + k + "/dipole" + str(polar)
                        criteria1.append(os.path.isfile(path + "/" + "xanes.dat"))
#                         criteria2.append(os.path.isfile(path + "/" + "xanes.dat"))
#                 criteria2 = list(set(criteria2))
#                 if not all(criteria2):
#                     eslst.remove(j)
            criteria1 = list(set(criteria1))
            if len(criteria1) == 1 and False in criteria1:
                gslst.remove(i)
        return gslst, eslst
    # cannot merge into screen_folder?
    def screen_folder2(gslst, eslst, absorber):
        for i in gslst:

            for j in eslst:
                criteria2 = []
                for k in absorber:
                    for polar in [1,2,3]:
                        path = folder + "/" + i + "/" + j + "/" + k + "/dipole" + str(polar)
#                         criteria1.append(os.path.isfile(path + "/" + "xanes.dat"))
                        criteria2.append(os.path.isfile(path + "/" + "xanes.dat"))
                criteria2 = list(set(criteria2))
                if not all(criteria2):
                    eslst.remove(j)

        return gslst, eslst

    gs = find_gs_folder()
    es = find_es_folder(gs)
    es, absorber = find_absorber(gs, es)
    gs_folder,es_folder = screen_folder(gs, es, absorber)
    gs_folder,es_folder = screen_folder2(gs, es, absorber)
    gs.sort()
    es.sort()
    absorber.sort()
    return gs, es, absorber

def XSAnaly(folder, gs, es, ab):
    spectra_objlst = []
    for i in gs:
        for j in es:
            path = folder + '/' + i + "/" +j
            spectra_objlst.append(XSplot(path, absorber=ab))

    # the last one is used as ref
    ref = spectra_objlst[-1]
    coss = []
    pearson = []
    spearman = []
    relArea = []
    alpha = []
    omega= []
    
    for i in spectra_objlst:

        comparePlots( i.spectra, ref.spectra, True, coss, pearson, spearman, relArea, omega, alpha )

    #reshape to 2D array
    coss = np.array(coss).reshape(len(gs),len(es))
    pearson = np.array(pearson).reshape(len(gs),len(es))
    spearman = np.array(spearman).reshape(len(gs),len(es))
    relArea = np.array(relArea).reshape(len(gs),len(es))

    # create dataframe
    coss_df=pd.DataFrame(coss, columns=es, index=gs)
    pearson_df=pd.DataFrame(pearson, columns=es, index=gs)
    spearman_df=pd.DataFrame(spearman, columns=es, index=gs)
    relArea_df=pd.DataFrame(relArea, columns=es, index=gs)
    
    # save dataframe to file
    coss_df.to_csv(folder+'_coss.csv')
    pearson_df.to_csv(folder+'_pearson.csv')
    spearman_df.to_csv(folder+'_spearman.csv')
    relArea_df.to_csv(folder+'_relArea.csv')

    # convert dataframe to inverse
    coss_df = coss_df.apply(lambda x: np.log10(1-x))
    pearson_df = pearson_df.apply(lambda x: np.log10(1-x))
    spearman_df = spearman_df.apply(lambda x: np.log10(1-x))
    relArea_df = relArea_df.apply(lambda x: np.log10(x))

    # save the new dataframe to file
    coss_df.to_csv(folder+'_coss_inv.csv')
    pearson_df.to_csv(folder+'_pearson_inv.csv')
    spearman_df.to_csv(folder+'_spearman_inv.csv')
    relArea_df.to_csv(folder+'_relArea_inv.csv')

def main():
    #print(sys.argv)
    if len(sys.argv) < 3:
        print('''
              Usage: code, folder
              code: anything starts with x stands for XSpectra; anything starts with o stands for OCEAN
              ''')
        exit()

    code = sys.argv[1]
    folder = sys.argv[2]
    
    if code.lower().startswith('x'):
        gs, es, ab = find_XS_folder(folder)
        XSAnaly(folder, gs, es, ab)
    elif code.lower().startswith('o'):
        print("OCEAN needs to be added")
    else:
        print("Code not supported")

if __name__ == '__main__':
    main()
