# ------------------------------------------------------------------------------
# Description:  Calculating wind speeds at different return periods
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-16
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
from helpers.pyExtras import getKeyList
from helpers.txtEditor import writeToTxt
import plotters.plot2D as plt

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

class response(object):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        """Scale time/freq. from model to full scale
        :param modelForce: obj w/ wind tunnel measurements
        :param RPeriod: str w/ return period of wind speed
        :param uH_fs, H_fs: flt w/ full scale building properties
        """
        self.modelForces = modelForces
        self.RPeriod= RPeriod
        self.H_fs   = H_fs
        self.uH_fs  = uH_fs
    
    def scaleTime(self):
        """Scale time/freq. from model to full scale
        """
        # Model properties
        self.modelForces.dT_ms    = 1 / self.modelForces.fq_sp_ms

        # Scaling factors
        self.lambda_u = self.uH_fs / self.modelForces.uH_ms
        self.lambda_g = self.H_fs / self.modelForces.H_ms
        self.lambda_fq= self.lambda_u / self.lambda_g
        self.lambda_t = 1 / self.lambda_fq

        # Scale quantities
        self.dT_fs    = self.lambda_t * self.modelForces.dT_ms
        self.fq_sp_fs   = self.lambda_fq * self.modelForces.fq_sp_ms

    def scaleForces(self):
        """Scale base forces from model to full scale
        """
        # Scaling factors
        self.lambda_F = self.lambda_u ** 2 * self.lambda_g ** 2
        self.lambda_M = self.lambda_u ** 2 * self.lambda_g ** 3

        # Scale floor forces
        self.modelForces.F_p_fs_D= self.modelForces.F_p_ms_D * self.lambda_F
        self.modelForces.F_p_fs_L= self.modelForces.F_p_ms_L * self.lambda_F
        
        # Scale base forces
        self.modelForces.BF_p_fs_D= self.modelForces.BF_p_ms_D * self.lambda_F
        self.modelForces.BF_p_fs_L= self.modelForces.BF_p_ms_L * self.lambda_F
        self.modelForces.BM_p_fs_D= self.modelForces.BM_p_ms_D * self.lambda_M
        self.modelForces.BM_p_fs_L= self.modelForces.BM_p_ms_L * self.lambda_M

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
        fq = np.fft.fftfreq(n, dT)[0:N]

        return S_p, fq

    def calcSpectralResponse(self, fq_p, S_p, fq_e, D):
        """Calculate the response spectrum
        """
        # Length of spectrum
        N = np.shape(S_p)[0]

        # Apply Dynamic amplification factor
        S_r = np.zeros(N)

        for i in range(0, N):
            eta_i   = fq_p[i]/fq_e
            H_i     = response.mechanicalAdmittance(self, eta_i, D)
            S_r[i]  = abs(H_i)**2 * S_p[i]

        return S_r

    def mechanicalAdmittance(self, eta, D):
        """Mechanical admittance of the 1-DOF System
        """
        # Dynamic amplification factor
        H_fq = 1 / np.sqrt((1-eta**2)**2 + (2*D*eta)**2)
        return H_fq

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

    def calcPeakFactor(self, fq_e, T):
        """Compute the peak factor
        """
        g_peak = np.sqrt(2 * np.log(fq_e * T)) + 0.5772 / np.sqrt((2 * np.log(fq_e * T)))
        return g_peak


class baseResponseForces(response):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        super().__init__(modelForces, RPeriod, uH_fs, H_fs)

    def calcResponse(self, fname, fq_e_D, fq_e_L, D):
        # Investigated wind speed
        writeToTxt(fname, "------------------------------")
        writeToTxt(fname, "u_H_mean:      " + '{:02.3f}'.format(self.uH_fs))
        writeToTxt(fname, "Return Period: " + self.RPeriod)
        writeToTxt(fname, "------------------------------")     

        # Base moment, drag direction
        writeToTxt(fname, "Base moment in drag direction [kNm]")
        baseResponseForces.calcLoadStats(self, fname, self.modelForces.BM_p_fs_D)
        baseResponseForces.calcPeakLoading(self, fname, self.modelForces.BM_p_fs_D, self.dT_fs, fq_e_D, D)
        writeToTxt(fname, "------------------------------")

        # Base moment, lift direction
        writeToTxt(fname, "Base moment in lift direction [kNm]")
        baseResponseForces.calcLoadStats(self, fname, self.modelForces.BM_p_fs_L)
        baseResponseForces.calcPeakLoading(self, fname, self.modelForces.BM_p_fs_L, self.dT_fs, fq_e_L, D)
        writeToTxt(fname, "------------------------------")

    def calcLoadStats(self, fname, F_p):
        # Calc statistics
        F_p_mean = np.mean(F_p)
        F_p_max  = np.max(F_p)
        F_p_min  = np.min(F_p)
        F_p_std  = np.std(F_p)

        writeToTxt(fname, "F_p_mean:      " + '{:02.3f}'.format(F_p_mean))
        writeToTxt(fname, "F_p_max:       " + '{:02.3f}'.format(F_p_max))
        writeToTxt(fname, "F_p_min:       " + '{:02.3f}'.format(F_p_min))
        writeToTxt(fname, "F_p_std:       " + '{:02.3f}'.format(F_p_std))

    def calcPeakLoading(self, fname, F_p, dT, fq_e, D):
        # Asses response with spectral analysis
        # --------------------  
        # Transform only the fluctuations "F_p_prime" to frequency domain
        F_p_mean  = np.mean(F_p)  
        F_p_prime = F_p - F_p_mean
        S_p, fq_p = response.transToSpectralDomain(self, F_p_prime, dT)

        # Apply mechanical admittance to the spectrum
        S_r = response.calcSpectralResponse(self, fq_p, S_p, fq_e, D)

        # # Setting up data to be plotted
        # plt.plot2D(fq_p, S_r, "f [Hz]", "Sr", "Spectrum", ["PSD"], xscale='log', yscale='log', savePlt=False, showPlt=True)

        # Perform the numerical integration
        F_r_std = response.numericalIntSpectrum(self, dT, S_r)

        # Estimate peak values
        g_peak = response.calcPeakFactor(self, 3600, fq_e)      # Peak Factor    
        F_r_max = F_p_mean + g_peak * F_r_std                   # Estimate max. response

        writeToTxt(fname, "Asses response with spectral analysis")
        writeToTxt(fname, "F_p_mean:        " + '{:02.3f}'.format(F_p_mean))
        writeToTxt(fname, "F_r_std:         " + '{:02.3f}'.format(F_r_std))
        writeToTxt(fname, "g_peak:          " + '{:02.3f}'.format(g_peak))
        writeToTxt(fname, "F_r_max:         " + '{:02.3f}'.format(F_r_max))

        # Comparison with the loading
        # --------------------  
        F_p_max = np.max(F_p)
        DLF_max = F_r_max / F_p_max

        writeToTxt(fname, "Comparison with loading")
        writeToTxt(fname, "F_p_max:         " + '{:02.3f}'.format(F_p_max))
        writeToTxt(fname, "DLF(F_max):      " + '{:02.3f}'.format(DLF_max))


class TipResponseDeflections(response):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        super().__init__(modelForces, RPeriod, uH_fs, H_fs)

class TipResponseAccelerations(response):
    def __init__(self, modelForces, RPeriod, uH_fs, H_fs):
        super().__init__(modelForces, RPeriod, uH_fs, H_fs)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     
