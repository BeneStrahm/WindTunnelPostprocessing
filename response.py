# ------------------------------------------------------------------------------
# Description:  
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-22
# Execution:    Import functions / collections (from folder.file import func)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  
import numpy as np
from scipy import integrate
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------
from helpers.txtEditor import writeToTxt
import plotters.plot2D as plt
from helpers.pyExtras import getKeyList
# ------------------------------------------------------------------------------
# Abbreviations
# ------------------------------------------------------------------------------
# p...  load
# r...  response
# ms... model scale
# fs... full scale 
# L...  lift
# D...  drag
# M...  moment
# F...  force
# H...  (at) height of building
# sp..  sample
# fq... frequency
# dn... direction

# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

def run_once(f):
        def wrapper(*args, **kwargs):
            if not wrapper.has_run:
                wrapper.has_run = True
                return f(*args, **kwargs)
        wrapper.has_run = False
        return wrapper

class amt():
    """Aerodynamic model theorie
    """
    def __init__(self):
        pass

    def transToSpectralDomain(self, F_p, dT):
        """Transform time series of forces into the spectral domain
        """
        # Length of time series
        n    = np.shape(F_p)[0]
        N    = n//2 

        # Get the Spectral Density, only positive half
        S_p = abs(np.fft.fft(F_p)[0:N])

        # Compute the power spectra density
        S_p = S_p ** 2 
        
        # According to "Boggs - Wind Loading ...[1991], (p.237)": S(fq)+ = 2 * S(fq)+/-
        S_p = 2 * S_p

        # Scaling Factor
        S_p = S_p / n

        # Scaling Factor
        S_p = S_p * dT 

        # Compute the frequencies, only positive half
        fq_p = np.fft.fftfreq(n, dT)[0:N]

        return S_p, fq_p

    def calcSpectralResponse(self, fq_p, S_p, fq_e, D, r="bm"):
        """Calculate the response spectrum
        :param r: response quantity (bm..base moment, a..acceleration)
        :type F_p: np.array[float]
        :param F_p: time series of aerodynamic forces
        :type F_p: np.array[float]
        """
        # Length of spectrum
        N = np.shape(S_p)[0]

        # Apply Dynamic amplification factor
        S_r = np.zeros(N)
        for i in range(0, N):
            eta_i   = fq_p[i]/fq_e
            if r == 'bm':
                H_i     = self.mechanicalAdmittance(eta_i, D)
                S_r[i]  = abs(H_i)**2 * S_p[i]
            elif r == 'a':
                Ha_i    = self.accelerationAdmittance(eta_i, fq_p[i], D)
                S_r[i]  = abs(Ha_i)** 2 * S_p[i]
            else:
                raise ('Invalid argument for desired response quantitiy')
        return S_r

    def mechanicalAdmittance(self, eta, D):
        """Mechanical admittance of the 1-DOF System
        """
        # Dynamic amplification factor
        H_fq = 1 / np.sqrt((1-eta**2)**2 + (2*D*eta)**2)
        return H_fq
    
    def accelerationAdmittance(self, eta, fq_p, D):
        """Acceleration admittance of the 1-DOF System
        """
        Ha_fq = (2 * np.pi * fq_p) ** 2 * self.mechanicalAdmittance(eta, D) 
        return Ha_fq

    def numericalIntSpectrum(self, dT, S_r):
        """Integrate the response spectrum
        """
        # Length of spectrum
        N = np.shape(S_r)[0]

        # Nyquist frequency
        fq_nyq = 1 / (2 * dT)

        # Sample spacing dfq
        dfq = fq_nyq / N 

        # Perform the numerical integration with Simpson rule
        F_ms = integrate.simps(S_r, dx=dfq) 

        # Return directly RMS
        F_rms = np.sqrt(F_ms)
        return F_rms

    def calcPeakFactor(self, nue, T_exp):
        """Compute the peak factor
        :param nue: effective cycling rate, freq. w/ most energy
        :type nue: float
        :param T_exp: exposure period, same basis as ref. wind speed
        :type T_exp: float
        """
        g_peak = np.sqrt(2 * np.log(nue * T_exp)) + 0.5772 / np.sqrt((2 * np.log(nue * T_exp)))
        return g_peak

    # @run_once
    def plotSpectrum(self, f, S, B, uH, F_ms, title, fname, vLines=None, vTexts=None, mode='reduced'):
        """Compute the peak factor
        :param nue: effective cycling rate, freq. w/ most energy
        :type nue: float
        :param T_exp: exposure period, same basis as ref. wind speed
        :type T_exp: float
        """
        # Function fit
        S_fit = savgol_filter(S, 101, 3)      # window size 51, polynomial order 3

        if mode == 'reduced':
            # Reduced frequency = f*B/U
            x = [f*B/uH]

            # Reduced PSD = f * S(f) / F_ms
            y = [f*S/F_ms, f*S_fit/F_ms]
            
            # Set up the labels
            xlabel = r'$f \cdot B \slash U_{H}$'
            ylabel = r'$f \cdot S(f) \slash \sigma^{2}$'

        elif mode == 'real':
            # Frequency
            x = [f]

            # PSD
            y = [S, S_fit]

            # Set up the labels
            xlabel = r'$f$'
            ylabel = r'$S(f)$'
            
        else:
            raise ('Invalid argument for spectral plot: Choose "reduced" or "real"')        

        # Crop to window with senseful data
        # xlim = [np.min(x[0])*10**1,np.max(x[0])*10**-1]
        ylim = [np.max(y[1])*10**-4,np.max(y[1])*10**1]

        legend = ["measured", "fitted"]
       
        style_dict = {"lines.linewidth":0.5, 'savefig.format':'svg'}
        
        plt.plot2D(x, y, xlabel, ylabel, title, legend, dir_fileName=fname,
                    vLines=vLines, vTexts=vTexts,  
                    xlim=[], ylim=ylim, xscale='log', yscale='log',
                    style_dict=style_dict, mpl='default', colorScheme='Monochrome', variation='color',
                    savePlt=True, savePkl=False, showPlt=False)

    # @run_once
    def plotLoadSpectum(self, windStats, buildProp, feModel, responseForces, fname, mode='reduced'):
        """Compute the peak factor
        :param nue: effective cycling rate, freq. w/ most energy
        :type nue: float
        :param T_exp: exposure period, same basis as ref. wind speed
        :type T_exp: float
        """
        f   = responseForces.fq_p
        S   = responseForces.S_p

        B   = buildProp.B
        uH  = buildProp.uH

        fq_e= feModel.fq_e

        uH_r= windStats.uH

        vLines, vTexts = [], []  

        for rPeriod in getKeyList(uH_r):
            if rPeriod in ["uH_050", "uH_002"]:
                vLines.append(fq_e * B / uH_r[rPeriod])
                vTexts.append(r'$R=$' + str(int(rPeriod[-3:])) + r'$ yr$')

        if buildProp.dn == 'D':
            title = "Load spectrum in drag direction"
        elif buildProp.dn == 'L':
            title = "Load spectrum in lift direction"

        F_ms=responseForces.F_p_std**2

        self.plotSpectrum(f, S, B, uH, F_ms, title, fname, vLines=vLines, vTexts=vTexts, mode=mode)
    
    # @run_once
    def plotResponseSpectrum(self, windStats, buildProp, feModel, responseForces, fname, mode='reduced'):
        """Compute the peak factor
        :param nue: effective cycling rate, freq. w/ most energy
        :type nue: float
        :param T_exp: exposure period, same basis as ref. wind speed
        :type T_exp: float
        """
        f   = responseForces.fq_p
        S   = responseForces.S_r

        B   = buildProp.B
        uH  = buildProp.uH

        fq_e= feModel.fq_e

        uH_r= windStats.uH

        if buildProp.dn == 'D':
            title = "Response spectrum in drag direction"
        elif buildProp.dn == 'L':
            title = "Response spectrum in lift direction"

        F_ms=responseForces.F_p_std**2

        self.plotSpectrum(f, S, B, uH, F_ms, title, fname, vLines=None, vTexts=None, mode=mode)
    

    def amtValidation(self):
        # --- Conventions --#
        # S_p:    : Spectrum of the loading  (p - loading)
        # S_r     : Spectrum of the response (r - response)
        # f       : Frequency of load / response spectra
        # H_f     : Mechanical Admittance of a certain frequency
        # D       : Damping ratio
        # omega_e : Circular eigenfrequency of the system
        
        # --- Input data ---#      
        # # Time settings
        dT = 1 / 2000
        eT = 20 
        nT = int(eT // dT + 1)
        T  = np.linspace(0, eT, nT)

        # Structural characteristis
        k       = (2 * np.pi)  ** 2                 #[kN/m]
        m       = 1                                 #[t]  
        omega_e = np.sqrt(k/m)
        f_e     = 1 / (2 * np.pi) * omega_e

        D       = 0.0125

        # Dynamic loading component
        F_G         = 1                                 #[kN] - self weight, can be neglected for dynamic amplification, see TUM: TM 3 - p.43
        F_0         = 1                                 #[kN]
        omega_p     = np.pi
        F_p         = F_G + F_0 * np.sin(omega_p * T)
        F_p_mean    = np.mean(F_p)   
        F_p_prime   = F_p - F_p_mean

        #######################
        # Analytical solution #
        #######################

        # Static deformation
        u_g    = F_G / k                                # [m]
        u_stat = F_0 / k                                # [m]

        # Frequency relation
        eta     = omega_p / omega_e   

        # Dynamic amplification factor
        V       = 1 / np.sqrt((1-eta**2)**2 + (2*D*eta)**2) 

        # Max Deformation
        u_r_max_ex = u_stat * V                            # [m]
        u_tot   = u_g + u_r_max_ex                         # [m]

        # Max Acceleration
        a_r_max_ex   = u_r_max_ex * omega_e ** 2                # [m/s]

        # Max response load
        F_r_max_ex = F_p_mean + V * np.max(F_p_prime)      # [kN]

        ######################
        # Aero Model Theorie #
        ######################
        # Transform only the fluctuations "F_p_prime" to frequency domain
        F_p_mean  = np.mean(F_p)  
        F_p_prime = F_p - F_p_mean

        # Transform time series of forces into the spectral domain
        S_p, fq_p = self.transToSpectralDomain(F_p_prime, dT)

        # Peak loading
        # Calculate the response spectrum
        S_r = self.calcSpectralResponse(fq_p, S_p, f_e, D, r='bm')

        # Integrate the response spectrum to get rms values
        F_r_rms = self.numericalIntSpectrum(dT, S_r)

        # Compute the peak factor, For Sinusoidal loading
        g_peak = np.sqrt(2)

        # Compute the peak loading
        F_r_max_amt = F_p_mean + g_peak * F_r_rms 

        # Peak accelerations
        # Calculate the response spectrum
        S_r = self.calcSpectralResponse(fq_p, S_p, f_e, D, r='a')

        # Integrate the response spectrum to get rms values
        a_r_std = self.numericalIntSpectrum(dT, S_r) / k
        
        # Compute the peak factor
        g_peak = np.sqrt(2)

        # Compute the peak acceleration
        a_r_max_amt = g_peak * a_r_std 

        print("Exact solution:")
        print("F_r_max: " + str(F_r_max_ex))
        print("u_r_max: " + str(u_r_max_ex))
        print("a_r_max: " + str(a_r_max_ex))

        print("Aerodynamic model theorie:")
        print("F_r_max: " + str(F_r_max_amt))
        print("a_r_max: " + str(a_r_max_amt))

class responseForces(amt):
    def __init__(self, F_p, dT, fq_e, D, nue, T_exp):
        """Inits the calculation acc. to the aerodynamic model theory
        :param F_p: time series of aerodynamic forces
        :type F_p: np.array[float]
        :param dT: time stepping of the time series [s]
        :type dT: float
        :param fq_e: eigenfrequency of the building [Hz]
        :type fq_e: float
        :param D: damping ratio [%]
        :type D: float
        :param nue: effective cycling rate, freq. w/ most energy
        :type nue: float
        :param T_exp: exposure period, same basis as ref. wind speed
        :type T_exp: float
        """
        super().__init__() 

        # Calc statistics
        self.F_p_mean = np.mean(F_p)
        self.F_p_max  = np.max(F_p)
        self.F_p_min  = np.min(F_p)
        self.F_p_std  = np.std(F_p)
        # Transform only the fluctuations "F_p_prime" to frequency domain
        self.F_p_prime = F_p - self.F_p_mean

        # Transform time series of forces into the spectral domain
        self.S_p, self.fq_p = super().transToSpectralDomain(self.F_p_prime, dT)

        # Integrate the load spectrum to get rms values
        self.F_p_std = super().numericalIntSpectrum(dT, self.S_p)

        # Calculate the response spectrum
        self.S_r = super().calcSpectralResponse(self.fq_p, self.S_p, fq_e, D, r='bm')

        # Integrate the response spectrum to get rms values
        self.F_r_std = super().numericalIntSpectrum(dT, self.S_r)

        # Compute the peak factor
        self.g_peak = super().calcPeakFactor(nue, T_exp)

        # Compute the peak loading
        self.F_r_max = self.F_p_mean + self.g_peak * self.F_r_std 
        self.F_r_min = self.F_p_mean - self.g_peak * self.F_r_std 

        # Comparison with the loading
        self.DLF_max = self.F_r_max / self.F_p_max

    def writeResultsToFile(self, fname, uH_fs, RPeriod):
        # Investigated wind speed
        writeToTxt(fname, "------------------------------")
        writeToTxt(fname, "u_H_mean:      " + '{:02.3f}'.format(uH_fs))
        writeToTxt(fname, "Return Period: " + RPeriod)
        writeToTxt(fname, "------------------------------")     

        # Load statistics
        writeToTxt(fname, "Aerodynamic loading")
        writeToTxt(fname, "F_p_mean:      " + '{:02.3f}'.format(self.F_p_mean))
        writeToTxt(fname, "F_p_max:       " + '{:02.3f}'.format(self.F_p_max))
        writeToTxt(fname, "F_p_std:       " + '{:02.3f}'.format(self.F_p_std))
        writeToTxt(fname, "------------------------------")  

        # Asses response with spectral analysis      
        writeToTxt(fname, "Asses response with spectral analysis")
        writeToTxt(fname, "F_r_std:         " + '{:02.3f}'.format(self.F_r_std))
        writeToTxt(fname, "g_peak:          " + '{:02.3f}'.format(self.g_peak))
        writeToTxt(fname, "F_r_max:         " + '{:02.3f}'.format(self.F_r_max))
        writeToTxt(fname, "------------------------------")  

        # Comparison with the loading
        writeToTxt(fname, "Comparison with loading")
        writeToTxt(fname, "DLF(F_max):      " + '{:02.3f}'.format(self.DLF_max))
        writeToTxt(fname, "------------------------------")  
    
class responseDeflection:
    """Class containing the full scale properties of the building
    :cvar coords: Cartesian coordinates of the node
    :vartype coords: list[float, float, float]
    """
    def __init__(self, feModel, responseForces, F_p_j):
        """Inits the class.
        :param coords: Cartesian coordinates of the node *([x], [x, y] or [x, y, z])*
        :type coords: list[float]
        """
        self.calcMeanDeflection(feModel, F_p_j)

        self.calcPeakDeflection(feModel, responseForces)

    def calcMeanDeflection(self, feModel, F_p_j):
        """Inits the class.
        :param coords: Cartesian coordinates of the node *([x], [x, y] or [x, y, z])*
        :type coords: list[float]
        """
        self.delta_tip_p_mean = feModel.calcStaticWindloadDeflection(F_p_j)

    def calcPeakDeflection(self, feModel, responseForces):
        """Inits the class.
        :param coords: Cartesian coordinates of the node *([x], [x, y] or [x, y, z])*
        :type coords: list[float]
        """
        # Calculate rms displacement
        self.delta_tip_r_std = feModel.v[0][0] / feModel.K_gen * responseForces.F_r_std

        # Get peak factor
        self.g_peak = responseForces.g_peak
        
        # Compute the peak loading
        self.delta_tip_r_max = self.delta_tip_p_mean + self.g_peak * self.delta_tip_r_std 
        self.delta_tip_r_min = self.delta_tip_p_mean - self.g_peak * self.delta_tip_r_std 
    
    def writeResultsToFile(self, fname, uH_fs, RPeriod):
        # Investigated wind speed
        writeToTxt(fname, "------------------------------")
        writeToTxt(fname, "u_H_mean:      " + '{:02.3f}'.format(uH_fs))
        writeToTxt(fname, "Return Period: " + RPeriod)
        writeToTxt(fname, "------------------------------")     

        # Load statistics
        writeToTxt(fname, "Aerodynamic loading")
        writeToTxt(fname, "delta_p_mean:    " + '{:02.3f}'.format(self.delta_tip_p_mean))
        writeToTxt(fname, "------------------------------")  

        # Asses response with spectral analysis      
        writeToTxt(fname, "Asses response with spectral analysis")
        writeToTxt(fname, "delta_r_std:     " + '{:02.3f}'.format(self.delta_tip_r_std ))
        writeToTxt(fname, "g_peak:          " + '{:02.3f}'.format(self.g_peak))
        writeToTxt(fname, "delta_r_max:     " + '{:02.3f}'.format(self.delta_tip_r_max))
        writeToTxt(fname, "------------------------------")  

class responseAccelerations(amt):
    def __init__(self, feModel, F_p, dT, fq_e, D, nue, T_exp):
        """Inits the calculation acc. to the aerodynamic model theory
        :param F_p: time series of aerodynamic forces
        :type F_p: np.array[float]
        :param dT: time stepping of the time series [s]
        :type dT: float
        :param fq_e: eigenfrequency of the building [Hz]
        :type fq_e: float
        :param D: damping ratio [%]
        :type D: float
        :param nue: effective cycling rate, freq. w/ most energy
        :type nue: float
        :param T_exp: exposure period, same basis as ref. wind speed
        :type T_exp: float
        """
        super().__init__() 

        # Transform only the fluctuations "F_p_prime" to frequency domain
        self.F_p_mean  = np.mean(F_p)  
        self.F_p_prime = F_p - self.F_p_mean

        # Transform time series of forces into the spectral domain
        self.S_p, self.fq_p = super().transToSpectralDomain(self.F_p_prime, dT)

        # Calculate the response spectrum
        self.S_r = super().calcSpectralResponse(self.fq_p, self.S_p, fq_e, D, r='a')

        # Integrate the response spectrum to get rms values
        self.a_r_std = super().numericalIntSpectrum(dT, self.S_r)
        self.a_r_std = self.a_r_std * feModel.v[0][0] / feModel.K_gen
        
        # Compute the peak factor
        self.g_peak = super().calcPeakFactor(nue, T_exp)

        # Compute the peak acceleration
        self.a_r_max = self.g_peak * self.a_r_std 

    def writeResultsToFile(self, fname, uH_fs, RPeriod):
        # Investigated wind speed
        writeToTxt(fname, "------------------------------")
        writeToTxt(fname, "u_H_mean:      " + '{:02.3f}'.format(uH_fs))
        writeToTxt(fname, "Return Period: " + RPeriod)
        writeToTxt(fname, "------------------------------")     

        # Asses response with spectral analysis      
        writeToTxt(fname, "Asses response with spectral analysis")
        writeToTxt(fname, "a_r_std:         " + '{:02.3f}'.format(self.a_r_std ))
        writeToTxt(fname, "g_peak:          " + '{:02.3f}'.format(self.g_peak))
        writeToTxt(fname, "a_r_max:         " + '{:02.3f}'.format(self.a_r_max))
        writeToTxt(fname, "------------------------------")  
    
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   

def validateAmt():
    aeroModelTheorie = amt()
    aeroModelTheorie.amtValidation()