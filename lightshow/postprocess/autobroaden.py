
"""Using known XAS exprimental spectrum and simulation to optimize the broadening pararmeters for other simulations"""


from scipy.stats import pearsonr as ps
from scipy.interpolate import interp1d
from scipy.stats import spearmanr
from scipy import optimize
from scipy.stats import norm
from scipy.stats import cauchy
from scipy.special import voigt_profile as voigt

import numpy as np


class Site_spectum:
    def __init__(
        self,
        xin,
        yin,
        VASP_out = [float, float, float],
        weight = 1.0,
        
    ):
        """
        The spectrum information of a each site
        Args:
            xin, yin: simulated spectrum for each site
            VASP_out: total energy of the final state with core hole, fermi energy of the final state, 
            and total energy of the initial state without core hole
            weight: the site weight in the crystal 
        """
        
        self.xin = xin
        self.yin = yin
        self.vasp_out = VASP_out
        self.Weight = weight
        
    @property
    def spectrum(self):
        return self.xin, self.yin
    
    @property
    def weight(self):
        return self.Weight
    
    @property
    def VASP_out(self):
        return self.vasp_out

    
class Broaden:
    def __init__(
        self,
        core_state = float,
        sigma = 1.0,
        gamma = 1.0,
        EXPshift = None
    ):
        """
        Args:
            core_state: Core state, e.g., 1s state energy for K-edge
            EXPshift: Empirical constant shift to match theory to experiment
        """
      
        self.core_state = core_state
        self.sigma = sigma
        self.gamma = gamma
        self.EXPshift = EXPshift

    
    @staticmethod
    def gauss_broaden(x,xin,yin,sigma):
        yout = np.zeros_like(x)
        
        num=len(xin)
        for i in range(num):

            yout=yout+yin[i]*norm.pdf(x, xin[i], sigma)/num   
        return yout 
    
    @staticmethod
    def lorentz_broaden(x,xin,yin,gamma):
        yout = np.zeros_like(x)
        
        num=len(xin)
        for i in range(num):

            yout=yout+yin[i]*cauchy.pdf(x, xin[i], gamma/2)/num   
        return yout 
   

    @staticmethod
    def voigt_broaden(x,xin,yin,sigma,gamma):
        yout = np.zeros_like(x)
        
        num=len(xin)
        for i in range(num):

            yout=yout+yin[i]*voigt(x-xin[i],sigma,gamma/2)/num   
        return yout 
    

    def crystal_spectrum_broaden(self, x, site_spectra, broad_method , sigma = None, gamma = None , shift = 0, cross_sec = True):
        
        yout = np.zeros_like(x)
        
        weight = 0
            
        for site_spectrum in site_spectra:
            xin, yin = site_spectrum.spectrum
        
            if self.EXPshift != None:
                xin = xin + self.EXPshift+site_spectrum.VASP_out[0]-site_spectrum.VASP_out[1]-site_spectrum.VASP_out[2]+self.core_state
            else:   
                xin = xin + shift
                
            if cross_sec:
                yin = yin * xin
                
            if broad_method == "gauss":
         
                yout=yout*weight + self.gauss_broaden(x,xin,yin,sigma)*site_spectrum.weight
            
            elif broad_method == "lorentz":

                yout=yout*weight + self.lorentz_broaden(x,xin,yin,gamma)*site_spectrum.weight
            
            elif broad_method == "voigt":
            
                yout=yout*weight + self.voigt_broaden(x,xin,yin,sigma,gamma)*site_spectrum.weight
            
            weight = weight + site_spectrum.weight            

        return yout/weight 
     
    
    
class Autobroaden():
    def __init__(
        self,
        broaden_paras: [float, float],
        CHlifetime : float,
        core_state : float,
        EXPshift = None
    ):

        """
        Args:
            broaden_paras: Broadening parameters, sigma of Gaussian and lorentz dividor 
            CHlifetime: Core hole life time, Zn K-edge: 1.67, Ti K-edge: 0.89, 
              taken from Campbell, J. and T. Papp, Widths of the atomic K–N7 levels. Atomic Data and Nuclear Data Tables, 2001. 77(1): p. 1-56.
            core_state: Core state, e.g., 1s state energy for K-edge
            EXPshift: Empirical constant shift to match theory to experiment
        """

        
        self.broaden_paras = broaden_paras
        self.CHlifetime = CHlifetime
        self.core_state = core_state
        self.EXPshift = EXPshift   

        
    @property
    def sigma(self):
        return self.broaden_paras[0]

    @property
    def lorentz_divider(self):
        return self.broaden_paras[1]
    
    
    @property
    def empirical_shift(self):
        return self.EXPshift 
    
  
    @staticmethod
    def deriv(x,y):
    
        yout = []

        for i in range(len(x)-1):
            yout.append((y[i+1]-y[i])/(x[i+1]-x[i]))

        yout.append(0)
        #yout = np.array(yout)
        
        return yout
        

    @staticmethod
    def voigt_broaden(x,xin,yin,sigma,ld,lt,ef):
        yout = np.zeros_like(x)
        
        num=len(xin)
        for i in range(num):
            gamma=lt/2
            if xin[i] > ef:
                gamma=gamma+(xin[i]-ef)/ld/2
            yout=yout+yin[i]*voigt(x-xin[i],sigma,gamma)/num   
        return yout 
    

    def crystal_spectrum_broaden(self, x, site_spectra, sigma, ld, shift , cross_sec = True):
        
        yout = np.zeros_like(x)
        
        weight = 0
          
        EXPshift = shift-site_spectra[0].VASP_out[0]+site_spectra[0].VASP_out[1]+site_spectra[0].VASP_out[2]-self.core_state

        for site_spectrum in site_spectra:
            xin, yin = site_spectrum.spectrum
                
            xin = xin + EXPshift+site_spectrum.VASP_out[0]-site_spectrum.VASP_out[1]-site_spectrum.VASP_out[2]+self.core_state
            
            if cross_sec:
                yin = yin * xin

            ef = EXPshift  + site_spectrum.VASP_out[0] -  site_spectrum.VASP_out[2] 
            yout=yout*weight + self.voigt_broaden(x,xin,yin,sigma,ld,self.CHlifetime,ef)*site_spectrum.weight
            
            weight = weight + site_spectrum.weight

        return yout/weight 
    
    
    def broaden_score(self,broaden_paras,x, site_spectra ,exp, dmu = True ,cross_sec = True, *args):
        
        spec = self.crystal_spectrum_broaden(x,site_spectra ,broaden_paras[0],broaden_paras[1], broaden_paras[2], cross_sec)
        
        if dmu:
            exp=self.deriv(x,exp) 
            spec=self.deriv(x,spec) 


        return -ps(spec,exp)[0]
    
    
    def broaden(self, x, site_spectra , shift = None, volume = 1.0, cross_sec = True):
        
        if shift == None:
            
            if self.EXPshift != None:
                
                shift = self.EXPshift+site_spectra[0].VASP_out[0]-site_spectra[0].VASP_out[1]-site_spectra[0].VASP_out[2]+self.core_state

            else:
                shift = 0
            
        return self.crystal_spectrum_broaden(x,site_spectra,self.broaden_paras[0], self.broaden_paras[1], shift, cross_sec)*volume


    
    
    def paras_optimize_scipy(
        self,
        x,
        site_spectra,
        exp,
        bounds = None,
        opt_shift = False,
        dmu = True,
        shift = 0,
        cross_sec = True
    ):
        """
        Using known experimental spectrum to optimize the broadening parameters 
        with simplicial homology global optimization method by scipy package
        
        """
        
        if opt_shift is False:
            
            xi = 10e-3
            
            bounds[2]=(shift-xi,shift+xi)
            
        result=optimize.shgo(self.broaden_score,bounds,args=(x,site_spectra,exp,dmu,cross_sec))

        self.broaden_paras = list(result.x) [0:2]
            
        if opt_shift is True:
            
            self.EXPshift = result.x[2] - site_spectra[0].VASP_out[0]+ site_spectra[0].VASP_out[1] + site_spectra[0].VASP_out[2] - self.core_state
            
        
        return result 