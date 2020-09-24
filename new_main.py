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

import aeroForces
import modelProp
import scaling
import fea
import response

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
    uH_f    = 46.10             # m/s
    H_f     = 128               # m
    f_e     = 46 / H_f          # Hz
    f_e     = 0.22706714
    D       = 0.02              # %
    I       = 477.924           # m4
    E       = 28900 * 10 ** 3   # kN/m2
    mue     = 30473 / H_f       # t/m

    # Load wind tunnel model properties, TPU Database files
    wtModelProp = modelProp.wtModelProp(fname)
    wtModelProp.loadTPUModelProp()

    # Initialize building model properties
    buildProp = modelProp.buildProp(H_f, E, I, I, mue, D, uH_f)
    
    # Load aerodynamic forces in model scale
    wtModelAeroForces = aeroForces.wtModelAeroForces()
    wtModelAeroForces.loadTPUModelForces(wtModelProp)

    # Calculate scaling factors
    scalingFactors = scaling.scalingFactors(wtModelProp, buildProp)

    # Scale wind tunnel model
    buildProp.scaleBuildProp(wtModelProp, scalingFactors)

    # Scale forces
    buildAeroForces = aeroForces.buildAeroForces(scalingFactors, wtModelAeroForces)

    # Create analysis model
    feModel = fea.feModel(buildProp, direction="L")

    # Compute eigenfrequencies
    feModel.getEigenfrequency()

    # Calc response forces
    responseForces = response.responseForces(buildAeroForces.BM_p_L, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)

    # Calc generalized quantities
    feModel.calcGeneralizedQuantitites()

    # Calc response deflections
    responseDeflections = response.responseDeflection(feModel, responseForces, buildAeroForces.F_p_L)

    # Calc response accelerations
    responseAccelerations = response.responseAccelerations(feModel, buildAeroForces.BM_p_L, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)

    print(responseAccelerations.a_r_max)

    print(responseDeflections.delta_tip_r_max * (feModel.fq_e * np.pi * 2) ** 2)







    # # Initialize class modelForce
    # modelForces = calcModelForces.modelForces(fname)

    # # Sum up base forces
    # modelForces.calcModelBaseForces()

    # # Sum up forces at each level
    # modelForces.calcModelFloorForces()

    # # Initialize class windStats
    # windStats = calcWindStats.windStats(uH_f)

    # # Calculate wind speeds at different return periods
    # windStats.calcWindStats()

    # # Calculate the base response forces for different wind speeds
    # for RPeriod in getKeyList(windStats.uH_fs):
    #     if "uH_fs_050" in RPeriod:
    #     # if True:
    #         baseResponseForces = calcResponse.baseResponseForces(modelForces, RPeriod, windStats.uH_fs[RPeriod], H_f)

    #         baseResponseForces.scaleTime()

    #         baseResponseForces.scaleForces()

    #         baseResponseForces.calcResponse(RPeriod + ".txt", f_e, f_e, D)

    # for RPeriod in getKeyList(windStats.uH_fs):
    #     if "uH_fs_050" in RPeriod:
    #     # if True:
    #         TipResponseDeflections = calcResponse.TipResponseDeflections(modelForces, RPeriod, windStats.uH_fs[RPeriod], H_f)

    #         TipResponseDeflections.scaleTime()

    #         TipResponseDeflections.scaleForces()

    #         TipResponseDeflections.calcResponse(RPeriod + ".txt", E, I, E, I, mue)

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
