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
# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

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
# fq...  frequency

# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

class amt():
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
                Ha_i     = self.accelerationAdmittance(eta_i, D, fq_p[i])
                S_r[i]  = abs(Ha_i)**2 * S_p[i]
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

    def amtValidation():
        pass


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

        # Transform only the fluctuations "F_p_prime" to frequency domain
        self.F_p_mean  = np.mean(F_p)  
        self.F_p_prime = F_p - self.F_p_mean

        # Transform time series of forces into the spectral domain
        self.S_p, self.fq_p = super().transToSpectralDomain(self.F_p_prime, dT)

        # Calculate the response spectrum
        self.S_r = super().calcSpectralResponse(self.fq_p, self.S_p, fq_e, D, r='bm')

        # Integrate the response spectrum to get rms values
        self.F_r_rms = super().numericalIntSpectrum(dT, self.S_r)

        # Compute the peak factor
        self.g_peak = super().calcPeakFactor(nue, T_exp)

        # Compute the peak loading
        self.F_r_max = self.F_p_mean + self.g_peak * self.F_r_rms 
        self.F_r_min = self.F_p_mean - self.g_peak * self.F_r_rms 

        # Comparison with the loading
        # --------------------  
        self.F_p_max = np.max(F_p)
        self.DLF_max = self.F_r_max / self.F_p_max
    
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
        self.delta_tip_r_std = feModel.v[0][0] / feModel.K_gen * responseForces.F_r_rms
        
        # Compute the peak loading
        self.delta_tip_r_max = self.delta_tip_p_mean + responseForces.g_peak * self.delta_tip_r_std 
        self.delta_tip_r_min = self.delta_tip_p_mean - responseForces.g_peak * self.delta_tip_r_std 

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
        self.a_r_rms = super().numericalIntSpectrum(dT, self.S_r)
        self.a_r_rms = self.a_r_rms * feModel.v[0][0] / feModel.K_gen
        
        # Compute the peak factor
        self.g_peak = super().calcPeakFactor(nue, T_exp)

        # Compute the peak acceleration
        self.a_r_max = self.g_peak * self.a_r_rms 
    
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   