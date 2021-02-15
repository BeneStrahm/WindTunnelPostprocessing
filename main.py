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
<<<<<<< HEAD
    fname = "C://Users//ac135564//GitHub//WindTunnelPostprocessing//T115_6//T115_6_000.mat"
=======
    ########
    # T115_6
    ########
    # fname = "C://Users//ac135564//GitHub//WindTunnelPostprocessing//T114_6//T114_6_000.mat"
    fname = "C://Users//bstra//GitHub//WindTunnelPostprocessing//T115_6//T115_6_000.mat"

    # File
    save_as= ('T115_6/results/SLS_design.txt')
>>>>>>> b97a13d8b6f35f1bb8a53369d5c218933b9d257f

    # Clean up results folder
    delFilesInFolder('T115_6/results')

    # Full scale building properties
<<<<<<< HEAD
    uH_f    = 37.59             # m/s
    H_f     = 160               # m
    B       = 32                # m
    D       = 0.02              # %
    I       = 477.924           # m4
    E       = 28900 * 10 ** 3   # kN/m2
    mue     = 30473 / H_f       # t/m
=======
    uH      = 38.96             # m/s       // Wind speed at z = H (50yr)
    H       = 160               # m         // Building height
    B       = 32                # m         // Building width
    nF      = 40                #           // Number of floors
    nM      = 5                 #           // Number of modules
    b       = 16                # m         // Core wall width 
    t       = 0.35              # m         // Core wall thickness   
    D       = 0.02              # %         // Damping
    I       = 956.191           # m4        // Starting value
    E       = 28900 * 10 ** 3   # kN/m2     // E-Modulus

    # Mass calculation (in t)
    M_DL_Floor  = 25600  
    M_DL_Col    = 401
    M_DL_IWall  = 3200

    # ########
    # # T114_6
    # ########
    # # fname = "C://Users//ac135564//GitHub//WindTunnelPostprocessing//T114_6//T114_6_000.mat"
    # fname = "C://Users//bstra//GitHub//WindTunnelPostprocessing//T114_6//T114_6_000.mat"

    # # File
    # save_as= ('T114_6/results/SLS_design.txt')

    # # Clean up results folder
    # delFilesInFolder('T114_6/results')

    # # Full scale building properties
    # uH      = 37.59             # m/s       // Wind speed at z = H (50yr)
    # H       = 128               # m         // Building height
    # B       = 32                # m         // Building width
    # nF      = 32                #           // Number of floors
    # nM      = 4                 #           // Number of modules
    # b       = 16                # m         // Core wall width 
    # t       = 0.175             # m         // Core wall thickness    
    # D       = 0.02              # %         // Damping
    # I       = 477.924           # m4        // Starting value
    # E       = 28900 * 10 ** 3   # kN/m2     // E-Modulus

    # # Mass calculation (in t)
    # M_DL_Floor  = 20480  
    # M_DL_Col    = 266
    # M_DL_IWall  = 2560

    #########################

    M_SDL_Floor = 1.0 * 0.100 * (B**2) * nF
    M_LL_Floor  = 0.2 * 0.250 * (B**2) * nF

>>>>>>> b97a13d8b6f35f1bb8a53369d5c218933b9d257f
    dns     = ['D', 'L']

    for dn in dns:    
        # Load wind tunnel model properties, TPU Database files
        wtModelProp = modelProp.wtModelProp(fname)
        wtModelProp.loadTPUModelProp()

        # Calculate wind stats at different return periods
        windStats = wind.windStats(uH)

        # Collect response quantities for plotting
        u_design = []
        M_r_max = []
        delta_r_max = []
        a_r_max = []

        for rPeriod in getKeyList(windStats.uH):
            if rPeriod in ["uH_050", "uH_002"]:
                # Initialize building model properties
                buildProp = modelProp.buildProp(H, B, nF, nM, dn, E, I, D, windStats.uH[rPeriod], \
                                                structSys='concreteCore', t=t, bCore=b)

                buildProp.calcBuildMass(M_DL_Floor, M_DL_Col, M_DL_IWall, M_SDL_Floor, M_LL_Floor)

                u_design.append(windStats.uH[rPeriod])

                
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
                responseForces.writeResultsToFile("T115_6/results/" + dn + "_forces.txt", windStats.uH[rPeriod], rPeriod)
                responseForces.plotLoadSpectum(windStats, buildProp, feModel, responseForces, "T115_6/results/" + dn + "_loadSpectrum", mode='reduced')
                responseForces.plotResponseSpectrum(windStats, buildProp, feModel, responseForces, "T115_6/results/" + dn + "_responseSpectrum", mode='real')

                # Calc response deflections
                responseDeflections = response.responseDeflection(feModel, responseForces, buildAeroForces.LF_p)
                responseDeflections.writeResultsToFile("T115_6/results/" + dn + "_deflections.txt", windStats.uH[rPeriod], rPeriod)
                delta_r_max.append(responseDeflections.delta_tip_r_max)

                # Calc response accelerations
                responseAccelerations = response.responseAccelerations(feModel, buildAeroForces.BM_p, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
                responseAccelerations.writeResultsToFile("T115_6/results/" + dn + "_accelerations.txt", windStats.uH[rPeriod], rPeriod)
                a_r_max.append(responseAccelerations.a_r_max)

        # plt.plot2D([u_design, u_design, u_design], [M_r_max/M_r_max[-1], delta_r_max/delta_r_max[-1] , a_r_max/a_r_max[-1]], 'Wind Speed', 'Normalized response',
        # 'Direction ' + dn, ['M_r_max', 'delta_r_max', 'a_r_max'], showPlt=True)

if __name__ == '__main__':
    main()
