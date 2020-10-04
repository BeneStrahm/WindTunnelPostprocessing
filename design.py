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
from helpers.txtEditor import writeToTxt

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------     

def main():
    # Get name of input file
    # fname = sys.argv [1]
    ########
    # T115_6
    ########
    # fname = "C://Users//ac135564//GitHub//WindTunnelPostprocessing//T114_6//T114_6_000.mat"
    fname = "C://Users//bstra//GitHub//WindTunnelPostprocessing//T115_6//T115_6_000.mat"

    # File
    save_as= ('T115_6/results/SLS_design.txt')

    # Clean up results folder
    delFilesInFolder('T115_6/results')

    # Full scale building properties
    uH      = 38.96             # m/s       // Wind speed at z = H (50yr)
    H       = 160               # m         // Building height
    B       = 32                # m         // Building width
    nF      = 40                #           // Number of floors
    nM      = 5                 #           // Number of modules
    b       = 16                # m         // Core wall thickness    
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
    # b       = 16                # m         // Core wall thickness    
    # D       = 0.02              # %         // Damping
    # I       = 477.924           # m4        // Starting value
    # E       = 28900 * 10 ** 3   # kN/m2     // E-Modulus

    # # Mass calculation (in t)
    # M_DL_Floor  = 20480  
    # M_DL_Col    = 266
    # M_DL_IWall  = 2560

    M_SDL_Floor = 1.0 * 0.100 * (B**2) * nF
    M_LL_Floor  = 0.2 * 0.250 * (B**2) * nF

    M_tot       = M_DL_Floor + M_DL_Col + M_DL_IWall + M_SDL_Floor + M_LL_Floor 

    mue     = (M_tot) / H       # t/m      # t/m       // Only from vertical loading/slabs 
    dns     = ['D', 'L']

    # Design criteria
    # Wall thickness
    t = np.arange(0.15, 0.85, 0.05)

    # Moment of Intertia   
    I = ((b+t)**4-(b-t)**4)/12

    # Stiff. red. per module
    stiff_red = 0.2

    # ULS
    # SLS
    u_lim   = H / 600         # u_r_max < u_lim
    a_lim   = 10 * 0.003        # a_r_rms < a_lim, ISO 10137 1-yr return limit for offices

    # Design for each direction
    for dn in dns:   
        
        writeToTxt(save_as, "------------------------------") 
        writeToTxt(save_as, "Direction " + dn)
        writeToTxt(save_as, "------------------------------") 

        # Load wind tunnel model properties, TPU Database files
        wtModelProp = modelProp.wtModelProp(fname)
        wtModelProp.loadTPUModelProp()

        # Calculate wind stats at different return periods
        windStats = wind.windStats(uH)
        
        # Design for deflections
        writeToTxt(save_as, "Design for deflections")
        # Set counter
        i = 0
        delta_tip_r_max = 10        
        while delta_tip_r_max > u_lim:
            t_iter = t[i]
            I_iter = I[i]
            rPeriod = "uH_050"  

            # Initialize building model properties
            buildProp = modelProp.buildProp(H, B, nF, nM, dn, E, I_iter,  mue, D, windStats.uH[rPeriod], \
                                            structSys='concreteCore', t=t_iter, bCore=b)
                
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
            feModel = fea.feModel(buildProp, stiff_red)

            # Compute eigenfrequencies
            feModel.getEigenfrequency()

            # Calc generalized quantities
            feModel.calcGeneralizedQuantitites()

            # Calc response forces
            responseForces = response.responseForces(buildAeroForces.BM_p, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)

            # Calc response deflections
            responseDeflections = response.responseDeflection(feModel, responseForces, buildAeroForces.LF_p)
            delta_tip_r_max = responseDeflections.delta_tip_r_max

            # Set counter for next iteration
            i = i + 1

        writeToTxt(save_as, "Iteration:           " + str(i))
        writeToTxt(save_as, "Wall thickness:      " + '{:02.3f}'.format(t[i]))
        writeToTxt(save_as, "Mom. of Intertia:    " + '{:02.3f}'.format(I[i]))
        writeToTxt(save_as, "Total mass:          " + '{:02.3f}'.format(feModel.MTot))
        writeToTxt(save_as, "Total mass:          " + '{:02.3f}'.format(np.sum(feModel.M)))
        writeToTxt(save_as, "u_r_max:             " + '{:02.3f}'.format(delta_tip_r_max))
        writeToTxt(save_as, "u_lim /u_r_max:      " + '{:02.3f}'.format(u_lim / delta_tip_r_max))
        writeToTxt(save_as, "--------------") 

        # Design for accelerations
        writeToTxt(save_as, "Design for accelerations")
        # Set counter
        i = 0
        a_r_rms = 10        
       
        while a_r_rms > a_lim:
            t_iter = t[i]
            I_iter = I[i]
            rPeriod = "uH_002" 
            
            # Initialize building model properties
            buildProp = modelProp.buildProp(H, B, nF, nM, dn, E, I_iter, mue, D, windStats.uH[rPeriod])
                
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
            feModel = fea.feModel(buildProp, t_iter)

            # Compute eigenfrequencies
            feModel.getEigenfrequency()

            # Calc generalized quantities
            feModel.calcGeneralizedQuantitites()

            # Calc response forces
            responseForces = response.responseForces(buildAeroForces.BM_p, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)

            # Calc response accelerations
            responseAccelerations = response.responseAccelerations(feModel, buildAeroForces.BM_p, buildProp.dT, feModel.fq_e, buildProp.D, feModel.fq_e, 360)
            a_r_rms = responseAccelerations.a_r_std

            # Set counter for next iteration
            i = i + 1

        writeToTxt(save_as, "Iteration:           " + str(i))
        writeToTxt(save_as, "Wall thickness:      " + '{:02.3f}'.format(t[i]))
        writeToTxt(save_as, "Mom. of Intertia:    " + '{:02.3f}'.format(I[i]))
        writeToTxt(save_as, "Total mass:          " + '{:02.3f}'.format(feModel.MTot))
        writeToTxt(save_as, "a_r_rms:             " + '{:02.3f}'.format(a_r_rms))
        writeToTxt(save_as, "a_lim /a_r_rms:      " + '{:02.3f}'.format(a_lim / a_r_rms))
        writeToTxt(save_as, "--------------") 

        # plt.plot2D([u_design, u_design, u_design], [M_r_max/M_r_max[-1], delta_r_max/delta_r_max[-1] , a_r_max/a_r_max[-1]], 'Wind Speed', 'Normalized response',
        # 'Direction ' + dn, ['M_r_max', 'delta_r_max', 'a_r_max'], showPlt=True)

if __name__ == '__main__':
    main()
