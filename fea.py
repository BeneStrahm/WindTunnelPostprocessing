# ------------------------------------------------------------------------------
# Description:  
# Author:       benedikt.strahm@ilek.uni-stuttgart.de
# Created:      2020-09-22
# Execution:    Import functions / collections (from folder.file import func)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------  

# ------------------------------------------------------------------------------
# Imported functions
# ------------------------------------------------------------------------------
import numpy as np
# For Static solver (feastruct)
from feastruct.pre.material import Material
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.solvers.feasolve import SolverSettings
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
# dn... direction
# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

class feModel:
    """Class containing finite element analyis of the building
    :cvar coords: Cartesian coordinates of the node
    :vartype coords: list[float, float, float]
    """
    def __init__(self, buildProp, stiff_red = 0):
        """Inits the class, setting up calculation model
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        """
        # Setting up calculation model
        # ------------
        # preprocessor
        # ---------
        # constants & lists
        self.L = buildProp.H                # length of the beam [m]
        self.n = buildProp.nz               # no of nodes [-]
        self.z = buildProp.z_lev            # coordinates [m]

        # everything starts with the analysis object
        self.analysis = FrameAnalysis2D()

        ### NODES
        # -------------
        # Coordinates for geometrie
        self.node_coords_geom = []
        for i in range(buildProp.nF+1):
            node = i * buildProp.hL
            self.node_coords_geom.append(node)
        
        # Sort from bottom to top of building, 
        # must be done for inserting loads at right order!
        self.node_coords_geom = sorted(self.node_coords_geom, reverse=True)

        self.node_coords_load = []
        for i in range(buildProp.nz):
            node = round(buildProp.z_lev[i], 1)
            self.node_coords_load.append(node)
        
        # Sort from bottom to top of building, 
        # must be done for inserting loads at right order!
        self.node_coords_load = sorted(self.node_coords_load, reverse=True)

        # Merge / sort / remove duplicate nodes
        self.node_coords = sorted(list(set(self.node_coords_geom + self.node_coords_load)), reverse=True)
            
        # Insert nodes 
        self.nodes = []
        for i in range(len(self.node_coords)):
            node = self.analysis.create_node(coords=[0, self.node_coords[i]])
            self.nodes.append(node)
        
        ## BEAMS
        # -------------
        self.beams = []

        # Counter for stiffness reduction
        j = -1

        # Determine coords. where reductions shall be made
        hModule = buildProp.H / buildProp.nM
        
        # Remark: Reverse node order here, going from bottom to tip
        for i in range(1, len(self.nodes)): #! n+1 (support, tip)
            # Get z-coord of current node
            z_i = self.nodes[-i].y

            # Raise counter if next module for reduction is reached
            if z_i % hModule == 0:
                j = j + 1

            # Recalculate dist. mass
            # For concrete walls: mue_walls = A_walls [m2] * 2.5 [t/m3]
            if buildProp.structSys == 'concreteCore':
                mue_walls = ((buildProp.bCore + buildProp.t)**2 - (buildProp.bCore - buildProp.t)**2) * 2.5 
                mue = buildProp.mue + (1 - stiff_red) ** j * mue_walls
            
            else:
                raise('No structural system specified')

            ## MATERIAL
            # -------------
            self.material = Material("Dummy", buildProp.E, 0.3, mue, colour='w')

            ## SECTION
            # -------------
            if buildProp.structSys == 'concreteCore':
                I = buildProp.I * (1 - stiff_red) ** j
            
            else:
                raise('No structural system specified')

            self.section = Section(area=1, ixx=I)

            # Add beams
            beam = self.analysis.create_element(
                   el_type='EB2-2D', nodes=[self.nodes[-i], self.nodes[-i-1]], material=self.material, section=self.section)
            self.beams.append(beam)
            
            # Print stiffness distribution for control
            # print('at z = ' + str(self.nodes[-i].y))
            # print('I = ' + str(beam.section.ixx))
            

        # boundary conditions are objects
        self.freedom_case = cases.FreedomCase()
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=0)
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=1)
        self.freedom_case.add_nodal_support(node=self.nodes[-1], val=0, dof=5)

    def getEigenfrequency(self):
        # ----------------
        # preprocessor
        # ----------------
        # an analysis case relates a support case to a load case
        analysis_case = cases.AnalysisCase(freedom_case=self.freedom_case, load_case=cases.LoadCase())

        # ----------------
        # frequency solver
        # ----------------
        settings = SolverSettings()
        settings.natural_frequency.time_info = True
        settings.natural_frequency.num_modes = 1

        solver = NaturalFrequency(
            analysis=self.analysis, analysis_cases=[analysis_case], solver_settings=settings)

        # Manual solver, see feastruct/solvers/naturalfrequency.py, in order
        # to extract mass/stiffness-matrix and eigenvectors      
        # assign the global degree of freedom numbers
        solver.assign_dofs()

        # Get the global stiffness / mass matrix
        (K, Kg) = solver.assemble_stiff_matrix()
        self.M = solver.assemble_mass_matrix()

        # apply the boundary conditions
        self.K_mod = solver.remove_constrained_dofs(K=K, analysis_case=analysis_case)
        self.M_mod = solver.remove_constrained_dofs(K=self.M, analysis_case=analysis_case)

        # Solve for the eigenvalues
        (self.fq_e, self.v) = solver.solve_eigenvalue(A=self.K_mod, M=self.M_mod, eigen_settings=settings.natural_frequency)

        # compute natural frequencies in Hz
        self.fq_e = np.sqrt(self.fq_e[0]) / 2 / np.pi

        # Normalize Eigenvector acc. to Boggs 1991, 4.3.1, p.109
        self.v = self.v / self.v[0] * self.L

        # Store stiffness matrix as np array
        self.K_mod = self.K_mod.toarray()
        self.M_mod = self.M_mod.toarray()

    def calcGeneralizedQuantitites(self):
        # Get generalized quantities
        self.K_gen = np.dot(np.dot(self.v.T, self.K_mod), self.v)[0][0]
        self.M_gen = np.dot(np.dot(self.v.T, self.M_mod), self.v)[0][0] 

        # To check, compute 
        # f_k = np.sqrt(K_gen/M_gen) / 2 / np.pi  d
        # print(f/f_k)
        # K_gen = 3 * E * I /(L^3) / 

    def calcStaticWindloadDeflection(self, F_p_j):
        # ------------
        # preprocessor
        # ---------

        # Adding loads
        self.load_case = cases.LoadCase()

        # Loop over all nodes
        for j in range(len(F_p_j)):
            # Search for node by z coordinate
            z_j = self.node_coords_load[j]
            for i in range(len(self.node_coords)):
                if z_j == self.node_coords[i]:
                    # Add loading to this node
                    F_p = np.mean(F_p_j[j])        #[in KN]
                    self.load_case.add_nodal_load(node=self.nodes[i], val=F_p , dof=0) 

        # an analysis case relates a support case to a load case
        analysis_case = cases.AnalysisCase(freedom_case=self.freedom_case, load_case=self.load_case)

        # ------
        # solver
        # ------

        # you can easily change the solver settings
        settings = SolverSettings()
        settings.linear_static.time_info = False

        # the linear static solver is an object and acts on the analysis object
        solver = LinearStatic(analysis=self.analysis, analysis_cases=[analysis_case], solver_settings=settings)
        solver.solve()

        # ----
        # post
        # ----
        # there are plenty of post processing options!
        # self.analysis.post.plot_geom(analysis_case=analysis_case)
        # self.analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=1e2)
        self.analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
        # self.analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
        # self.analysis.post.plot_reactions(analysis_case=analysis_case)
        
        # Support reactions, to check bending moment for validation
        for support in analysis_case.freedom_case.items:
            if support.dof in [5]:
                reaction = support.get_reaction(analysis_case=analysis_case)
        
        print('BASE FORCE')
        print(reaction)

        # read out deformation at top 
        delta_tip = self.nodes[0].get_displacements(analysis_case)[0]
       
        return delta_tip
    

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   