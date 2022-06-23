import numpy as np
from pathlib import Path
## pasers for XS, OCEAN and EXCITING
class XSplot():
    def __init__(self, folder, absorber=['0'], polar=['1', '2', '3']):
        self.path = Path(folder)
        self.absorber = absorber
        self.polar = polar

    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            for polar in self.polar:
                dipole = "dipole" + polar
                path = self.path / iab / dipole
                #print(path)
                if Spectra is None:
#                     Spectra = np.loadtxt(path + "/xanes.dat", comments='#')
                    Spectra = np.loadtxt( path / "xanes.dat", skiprows=4, usecols=(0,1)  )
                else:
                    Spectra[:,1] += np.loadtxt( path / "xanes.dat", skiprows=4, usecols=(1)  )
#                     Spectra[:,1] += np.loadtxt(path + "/xanes.dat", comments='#') # did not take average
        return Spectra

    def exists(self):
        out = []
        for iab in self.absorber:
            for polar in self.polar:
                dipole = "dipole" + polar
                path = self.path  /  iab  /  dipole
                if not Path(path / "xanes.dat").exists():
                    out.append(False)
                    print(f'{Path(path / "xanes.dat")} not exists')
                else:
                    out.append(True)
        if all(out):
            return True
        else:
            return False
                    

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


class XSplot_rescale(XSplot):
    def __init__(self, folder, volume, egrid, e1s, ehomo, absorber=['0'], polar=['1', '2', '3']):
        super(XSplot_rescale, self).__init__(folder, absorber, polar)
        self.egrid = np.linspace(egrid[0], egrid[1], egrid[2])
        self.e1s = e1s
        self.homo = ehomo
        self.omega = self.egrid + self.homo - self.e1s
        self.omega_ry = self.omega/13.6056980659/2 # 2 is from ry to hartree used in OCEAN & EXCITING
        self.volume = volume
    @property
    def y(self):
        return self.spectra[:,1] /self.omega_ry*137.04/self.volume
    @property
    def ny(self):
        return self.y
    @property
    def scaspectra(self):
        return np.stack((self.x,self.y)).T

class OCEANplot():
    def __init__(self, folder, absorber=['0'],polar = ['1', '2', '3']):
        self.path = Path(folder)
        self.absorber = absorber
        self.polar = polar

    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            for polar in self.polar:
                element = "Ti" # hard coded
                p = int(polar)
                site = int(iab) + 1
                path = self.path / "absspct_{:s}.{:04d}_1s_{:02d}".format( element, site, p)
                #print(path)
                if Spectra is None:
#                     Spectra = np.loadtxt(path + "/xanes.dat", comments='#')
                    Spectra = np.loadtxt( path, skiprows=2, usecols=(0,2)  )
                else:
                    Spectra[:,1] += np.loadtxt( path, skiprows=2, usecols=(2)  )
#                     Spectra[:,1] += np.loadtxt(path + "/xanes.dat", comments='#') # did not take average
        return Spectra

    def exists(self):
        out = []
        element = "Ti"
        for iab in self.absorber:
            for polar in self.polar:
                p = int(polar)
                site = int(iab) + 1
                path = self.path / "absspct_{:s}.{:04d}_1s_{:02d}".format( element, site, p)

                if not path.exists():
                    out.append(False)
                    print(f'{path} not exists')
                else:
                    out.append(True)
        if all(out):
            return True
        else:
            return False

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
    def __init__(self, folder, absorber=['0'], polar = ['11', '22', '33'], ip=False):
        self.path = Path(folder)
        self.absorber = absorber
        self.polar = polar
        self.ip = ip
    @property
    def spectra(self):
        Spectra = None
        for iab in self.absorber:
            for polar in self.polar:
                if not self.ip:
                    path = self.path  / iab / f"EPSILON/EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC{polar}.OUT" #.format( element, site, p)
                else:
                    path = self.path  / iab / f"EPSILON/EPSILON_BSE-IP_SCR-full_OC{polar}.OUT"

                if Spectra is None:
    #                     Spectra = np.loadtxt(path + "/xanes.dat", comments='#')
                    Spectra = np.loadtxt( path, skiprows=18, usecols=(0,2)  )
                else:
                    Spectra[:,1] += np.loadtxt( path, skiprows=18, usecols=(2)  )
        return Spectra

    def exists(self):
        out = []
        for iab in self.absorber:
            for polar in self.polar:
                dipole = "dipole" + polar
                if not self.ip:
                    path = self.path  /  iab  /  f"EPSILON/EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC{polar}.OUT"
                else:
                    path = self.path  / iab / f"EPSILON/EPSILON_BSE-IP_SCR-full_OC{polar}.OUT"
                if not path.exists():
                    out.append(False)
#                     print(f'{Path(path / "xanes.dat")} not exists')
                else:
                    out.append(True)
        if all(out):
            return True
        else:
            return False


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

