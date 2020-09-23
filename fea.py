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

# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------  

class feModel:
    """Class containing finite element analyis of the building
    :cvar coords: Cartesian coordinates of the node
    :vartype coords: list[float, float, float]
    """
    def __init__(self, buildProp, direction="D"):
        """Inits the class, setting up calculation model
        :param buildProp: full scale properties of the building
        :type buildProp: :class:`~modelProp.buildProp`
        :param direction: Direction of model [D,L]
        :type direction: :string "D" / "L"
        """
        self.direction = direction
        # Setting up calculation model
        # ------------
        # preprocessor
        # ---------
        # constants & lists
        self.L = buildProp.H                # length of the beam [m]
        self.n = buildProp.nz               # no of nodes [-]
        self.z = buildProp.z_lev            # coordinates [m]
        self.z = np.append(self.L, self.z)  # coordinates [m], append top of building
        self.z = np.append(self.z, 0)       # coordinates [m], append support node

        # everything starts with the analysis object
        self.analysis = FrameAnalysis2D()

        # materials and sections are objects
        self.material = Material("Dummy", buildProp.E, 0.3, buildProp.mue, colour='w')
        if self.direction == "D":
            self.section = Section(area=1, ixx=buildProp.Idd)
        elif self.direction == "L":
            self.section = Section(area=1, ixx=buildProp.Ill)
        else:
            raise Exception("Invalid direction specified") 

        # nodes are objects
        self.nodes = []
        for i in range(0,self.n+2): #! n+2 (support, tip)
            node = self.analysis.create_node(coords=[0, self.z[i]])
            self.nodes.append(node)

        # and so are beams!
        self.beams = []
        for i in range(0,self.n+1): #! n+1 (support, tip)
            beam = self.analysis.create_element(
                el_type='EB2-2D', nodes=[self.nodes[i], self.nodes[i+1]], material=self.material, section=self.section)
            self.beams.append(beam)

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
        M = solver.assemble_mass_matrix()

        # apply the boundary conditions
        self.K_mod = solver.remove_constrained_dofs(K=K, analysis_case=analysis_case)
        self.M_mod = solver.remove_constrained_dofs(K=M, analysis_case=analysis_case)

        # Solve for the eigenvalues
        (self.f, self.v) = solver.solve_eigenvalue(A=self.K_mod, M=self.M_mod, eigen_settings=settings.natural_frequency)

        # compute natural frequencies in Hz
        self.f = np.sqrt(self.f[0]) / 2 / np.pi

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
        

    def calcStaticMeanWindloadDeflection(self, buildAeroForces):
        # ------------
        # preprocessor
        # ---------

        # Adding loads
        self.load_case = cases.LoadCase()

        for i in range(self.n):
            if self.direction == "D":
                F_p = np.mean(buildAeroForces.F_p_D[i])        #[in KN]
            elif self.direction == "L":
                F_p = np.mean(buildAeroForces.F_p_L[i])        #[in KN]
            else:
                raise Exception("Invalid direction specified") 
            self.load_case.add_nodal_load(node=self.nodes[i+1], val=F_p , dof=0) # i+1 (support, tip)

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
        # analysis.post.plot_geom(analysis_case=analysis_case)
        # analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=1e2)
        # analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
        # analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
        # analysis.post.plot_reactions(analysis_case=analysis_case)
        
        # Support reactions, to check bending moment for validation
        for support in analysis_case.freedom_case.items:
            if support.dof in [5]:
                reaction = support.get_reaction(analysis_case=analysis_case)

        # read out deformation at top 
        self.delta_p_mean = self.nodes[0].get_displacements(analysis_case)[0]
       
    

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------   