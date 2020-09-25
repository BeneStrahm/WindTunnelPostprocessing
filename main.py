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
from helpers.filemanager import delFilesInFolder

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     

def main():
    # Get name of input file
    # fname = sys.argv [1]
    fname = "C://Users//ac135564//GitHub//WindTunnelPostprocessing-1//T114_4_000_1.mat"

    # Clean up results folder
    delFilesInFolder('results')

    # Full scale building properties
    uH_f    = 46.10             # m/s
    H_f     = 128               # m
    D       = 0.02              # %
    I       = 477.924           # m4
    E       = 28900 * 10 ** 3   # kN/m2
    mue     = 30473 / H_f       # t/m
    dns     = ['D', 'L']

    for dn in dns:    
        # Load wind tunnel model properties, TPU Database files
        wtModelProp = modelProp.wtModelProp(fname)
        wtModelProp.loadTPUModelProp()

        # Calculate wind stats at different return periods
        windStats = wind.windStats(uH_f)

        # Collect response quantities for plotting
        u_design = []
        M_r_max = []
        delta_r_max = []
        a_r_max = []

        for RPeriod in getKeyList(windStats.uH):
            if "uH_050" in RPeriod:
                # Initialize building model properties
                buildProp = modelProp.buildProp(H_f, dn, E, I, mue, D, windStats.uH[RPeriod])
                u_design.append(windStats.uH[RPeriod])
                
                # Load aerodynamic forces in model scale
                wtModelAeroForces = aeroForces.wtModelAeroForces()
                wtModelAeroForces.loadTPUModelForces(wtModelProp, buildProp)

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

                # # Calc response forces
                responseForces = response.responseForces(buildAeroForces.BM_p, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
                responseForces.writeResultsToFile("results/" + dn + "_forces.txt", windStats.uH[RPeriod], RPeriod)
                M_r_max.append(responseForces.F_r_max)

                # Calc response deflections
                responseDeflections = response.responseDeflection(feModel, responseForces, buildAeroForces.LF_p)
                responseDeflections.writeResultsToFile("results/" + dn + "_deflections.txt", windStats.uH[RPeriod], RPeriod)
                delta_r_max.append(responseDeflections.delta_tip_r_max)

                # Calc response accelerations
                responseAccelerations = response.responseAccelerations(feModel, buildAeroForces.BM_p, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
                responseAccelerations.writeResultsToFile("results/" + dn + "_accelerations.txt", windStats.uH[RPeriod], RPeriod)
                a_r_max.append(responseAccelerations.a_r_max)

        # plt.plot2D([u_design, u_design, u_design], [M_r_max/M_r_max[-1], delta_r_max/delta_r_max[-1] , a_r_max/a_r_max[-1]], 'Wind Speed', 'Normalized response',
        # 'Direction ' + dn, ['M_r_max', 'delta_r_max', 'a_r_max'], showPlt=True)

if __name__ == '__main__':
    main()
