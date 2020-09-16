# ------------------------------------------------------------------------------
# Project:      Wind tunnel postprocessing
# Description:  Main file
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-15
# Execution:    Drag & drop file on script
#               https://mindlesstechnology.wordpress.com/2008/03/29/make-python-scripts-droppable-in-windows/
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  

import sys           
import numpy as np                                             

# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------

import calcModelForces 
import calcWindStats
import calcResponse

from helpers.pyExtras import getKeyList
# import calcForcesFromWindtunnel as FFWT

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     

def main():
    # Get name of input file
    fname = sys.argv [1]

    # Full scale building properties
    uH_f    = 41.13
    H_f     = 160.0 
    f_e     = 46 / H_f
    omega_e = 2 * np.pi * f_e
    D       = 0.02

    modelForces = calcModelForces.modelForces(fname)

    modelForces.loadWindtunnelData()

    modelForces.calcModelBaseForces()

    modelForces.calcModelFloorForces()

    windStats = calcWindStats.windStats(uH_f)

    windStats.calcWindStats()

    for key, value in getKeyList(windStats.uH_fs):
        print(key)

    baseResponse = calcResponse.response(modelForces, uH_f, H_f)

    baseResponse.scaleTime()

    baseResponse.scaleForces()

    # # Calculate Forces from wind tunnel results
    # Fpm_D, Fpm_L, Mpm_D, Mpm_L = FFWT.calcModelBaseForces(filename)

    # # Scale time / forces from wind tunnel results
    # Fpf_D, Fpf_L, Mpf_D, Mpf_L = FFWT.scaleForces(filename, Fpm_D, Fpm_L, Mpm_D, Mpm_L, uH_f, H_f)
    # dT_f, nT, f_f = FFWT.scaleTime(filename, uH_f, H_f)

    # # Calculate full scale structural response force with the aerodynamic model theory
    # # --- Calculation --# 
    # S_p, S_r, f = AMT.SpectralAnalyis(Fpf_D, dT_f, omega_e, D)

    # #---- Plotting ----#
    # if False:
    #     # Plotting the load spectrum
    #     yLabel  = r'\$ S_p(f) \$'
    #     title   = "Load spectrum"
    #     AMT.plotSpectra(f, S_p, yLabel, title)


    #     # Plotting the response spectrum
    #     yLabel  = r'\$ S_r(f) \$'
    #     title   = "Response spectrum"
    #     AMT.plotSpectra(f, S_r, yLabel, title)

    # # wait for enter, otherwise we'll just close on exit
    input()

if __name__ == '__main__':
    main()
