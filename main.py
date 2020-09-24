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
import wind

import plotters.plot2D as plt
from helpers.pyExtras import getKeyList

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
    D       = 0.02              # %
    I       = 477.924           # m4
    E       = 28900 * 10 ** 3   # kN/m2
    mue     = 30473 / H_f       # t/m

    # Load wind tunnel model properties, TPU Database files
    wtModelProp = modelProp.wtModelProp(fname)
    wtModelProp.loadTPUModelProp()

    # Calculate wind stats at different return periods
    windStats = wind.windStats(uH_f)

    for RPeriod in getKeyList(windStats.uH):
        # if "uH_050" in RPeriod:
            # Initialize building model properties
            buildProp = modelProp.buildProp(H_f, E, I, mue, D, windStats.uH[RPeriod])
            
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
            feModel = fea.feModel(buildProp)

            # Compute eigenfrequencies
            feModel.getEigenfrequency()

            # Calc generalized quantities
            feModel.calcGeneralizedQuantitites()

            # Calc response forces
            responseForcesD = response.responseForces(buildAeroForces.BM_p_D, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
            responseForcesD.writeResultsToFile("results/Drag_forces.txt", windStats.uH[RPeriod], RPeriod)

            responseForcesL = response.responseForces(buildAeroForces.BM_p_L, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
            responseForcesL.writeResultsToFile("results/Lift_forces.txt", windStats.uH[RPeriod], RPeriod)

            # Calc response deflections
            responseDeflectionsD = response.responseDeflection(feModel, responseForcesD, buildAeroForces.F_p_D)
            responseDeflectionsL = response.responseDeflection(feModel, responseForcesL, buildAeroForces.F_p_L)

            # Calc response accelerations
            responseAccelerationsD = response.responseAccelerations(feModel, buildAeroForces.BM_p_D, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
            responseAccelerationsL = response.responseAccelerations(feModel, buildAeroForces.BM_p_L, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)

if __name__ == '__main__':
    main()
