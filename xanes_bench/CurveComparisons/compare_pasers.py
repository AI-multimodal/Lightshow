import numpy as np
## pasers for XS, OCEAN and EXCITING
class XSplot():
    def __init__(self, folder, absorber=['0'],polar=['1', '2', '3']):
        self.path = folder
        self.absorber = absorber
        self.polar=polar

    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            for polar in self.polar:
                path = self.path + "/" + iab + "/" + "dipole" + polar
                #print(path)
                if Spectra is None:
#                     Spectra = np.loadtxt(path + "/xanes.dat", comments='#')
                    Spectra = np.loadtxt( path + "/xanes.dat", skiprows=4, usecols=(0,1)  )
                else:
                    Spectra[:,1] += np.loadtxt( path + "/xanes.dat", skiprows=4, usecols=(1)  )
#                     Spectra[:,1] += np.loadtxt(path + "/xanes.dat", comments='#') # did not take average
        return Spectra

    @property
    def x(self):
        return self.spectra[:,0]
    @property
    def y(self):
        return self.spectra[:,1]
    @property
    def nspectra(self):
        x = self.x
        y = self.y/np.max(self.y)
        nplot = np.vstack((x,y))
        return nplot.T
    @property
    def nx(self):
        return self.nspectra[:,0]
    @property
    def ny(self):
        return self.nspectra[:,1]


class OCEANplot():
    def __init__(self, folder, absorber=['0'],polar = ['1', '2', '3']):
        self.path = folder
        self.absorber = absorber
        self.polar = polar

    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            for polar in self.polar:
                element = "Ti"
                p = int(polar)
                site = int(iab) + 1
                path = self.path + "/" + "absspct_{:s}.{:04d}_1s_{:02d}".format( element, site, p)
                #print(path)
                if Spectra is None:
#                     Spectra = np.loadtxt(path + "/xanes.dat", comments='#')
                    Spectra = np.loadtxt( path, skiprows=2, usecols=(0,2)  )
                else:
                    Spectra[:,1] += np.loadtxt( path, skiprows=2, usecols=(2)  )
#                     Spectra[:,1] += np.loadtxt(path + "/xanes.dat", comments='#') # did not take average
        return Spectra
    @property
    def x(self):
        return self.spectra[:,0]
    @property
    def y(self):
        return self.spectra[:,1]
    @property
    def nspectra(self):
        x = self.x
        y = self.y/np.max(self.y)
        nplot = np.vstack((x,y))
        return nplot.T
    @property
    def nx(self):
        return self.nspectra[:,0]
    @property
    def ny(self):
        return self.nspectra[:,1]

class Excitingplot():
    def __init__(self, folder, polar = ['11', '22', '33']):
        self.path = folder
#         self.absorber = absorber
        self.polar = polar
        
    @property    
    def spectra(self):
        Spectra = None
#         for iab in self.absorber:
        for polar in self.polar:
#             element = "Ti"
#             p = int(polar)
#             site = int(iab) + 1
            path = self.path  + f"/EPSILON_BSE-IP_SCR-full_OC{polar}.OUT" #.format( element, site, p) 
            #print(path)
            if Spectra is None:
#                     Spectra = np.loadtxt(path + "/xanes.dat", comments='#')
                Spectra = np.loadtxt( path, skiprows=18, usecols=(0,2)  )
            else:
                Spectra[:,1] += np.loadtxt( path, skiprows=18, usecols=(2)  )
#                     Spectra[:,1] += np.loadtxt(path + "/xanes.dat", comments='#') # did not take average
        return Spectra
    @property
    def x(self):
        return self.spectra[:,0] - self.spectra[0,0] - 15
    @property
    def y(self):
        return self.spectra[:,1]
    @property
    def nspectra(self):
        x = self.x
        y = self.y/np.max(self.y)
        nplot = np.vstack((x,y))
        return nplot.T
    @property
    def nx(self):
        return self.nspectra[:,0]
    @property
    def ny(self):
        return self.nspectra[:,1]