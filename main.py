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
    # fname = sys.argv [1]
    fname = "C://Users//ac135564//GitHub//WindTunnelPostprocessing-1//T115_4_000_1.mat"

    # Full scale building properties
    uH_f    = 46.10
    H_f     = 128
    f_e     = 46 / H_f
    D       = 0.02

    # Initialize class modelForce
    modelForces = calcModelForces.modelForces(fname)

    # Load wind tunnel data
    modelForces.loadWindtunnelData()

    # Sum up base forces
    modelForces.calcModelBaseForces()

    # Sum up forces at each level
    modelForces.calcModelFloorForces()

    # Initialize class windStats
    windStats = calcWindStats.windStats(uH_f)

    # Calculate wind speeds at different return periods
    windStats.calcWindStats()

    # Calculate the base response forces for different wind speeds
    for RPeriod in getKeyList(windStats.uH_fs):
        baseResponseForces = calcResponse.baseResponseForces(modelForces, RPeriod, windStats.uH_fs[RPeriod], H_f)

        baseResponseForces.scaleTime()

        baseResponseForces.scaleForces()

        baseResponseForces.calcResponse(RPeriod + ".txt", f_e, f_e, D)

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
    # input()

if __name__ == '__main__':
    main()
