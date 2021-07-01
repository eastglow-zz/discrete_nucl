#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100 #Number of elements in the x-direction
  ny = 100 #Number of elements in the y-direction
  xmax = 100 #X-direction domain size, for this file it is in nm
  ymax = 100 #Y-direction domain size, for this file it is in nm
[]

[GlobalParams] #Parameters used in multiple input blocks
  fa_name  = fb
  fb_name  = fv
  ca       = cvB
  cb       = cvV
[]

[Variables]
  # vacancy concentration
  [./cv]
  [../]
[]

[AuxVariables]
  [./f_dens]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./psi_map]
    order = CONSTANT
    family = MONOMIAL
  [../]
  #Identifies surface regions for heterogeneous nucleation
  [./surf]
  [../]

  # order parameter
  [./psi]
  [../]

  # bulk phase concentration (matrix)
  [./cvB]
  [../]

  # hydrogen phase concentration (delta phase)
  [./cvV]
  [../]

  [./cv_avg]
  [../]
[]

[ICs]
  [./cv_IC]
    type = ConstantIC
    variable = cv
    value = 2.021e-7
  [../]
  [./psi_IC]
    type = ConstantIC
    variable = psi
    value = 0
  [../]
  [./cvV_IC]
    type = ConstantIC
    variable = cvV
    value = 1.0
  [../]
  [./cvB_IC]
    type = ConstantIC
    variable = cvB
    value = 2.021e-7
  [../]

  # [./psi_circle_IC]
  #   type = SmoothCircleIC
  #   variable = psi
  #   radius = 2
  #   int_width = 1
  #   profile = TANH
  #   invalue = 1
  #   outvalue = 0
  #   x1 = 50
  #   y1 = 50
  # [../]

  [./surf_IC]
    type = BoundingBoxIC
    variable = surf
    x1 = -1
    x2 = 2
    y1 = -1
    y2 = 101
    inside = 1
    outside = 0
    int_width = 0.3
  [../]
[]

[Materials]
  #######Parameters defining the material behavior
  [./Conditions] #Sets the temperature, length scale, and time scale of the simulation
    type = GenericConstantMaterial
    prop_names =  'T   length_scale time_scale'
    prop_values = '300 1e-9         1e-6      '
  [../]
  [./Material_properties] #Sets the free energy properties of the material
    type = GenericConstantMaterial
    prop_names =  'Em    Ef   Ab    Bb eps  Av   cVeq D0'
    prop_values = '0.055 0.52 270.5 32 0.01 1000 1    8.37707e-7 '
  [../]
  [./Model_parameters]
    type = GenericConstantMaterial
    prop_names =  'L0   delta'
    prop_values = '0.1  0.6928'
  [../]
  [./Constants]
    type = GenericConstantMaterial
    prop_names =  'eVpJ           kb             pi'
    prop_values = '6.24150934e+18 8.6173303e-5   3.141592'
  [../]

  [./Va]
    type = ParsedMaterial
    f_name = Va
    # material_property_names = 'length_scale'
    # function = '5.586e24*length_scale^3'
    function = '2.24e-2'  #nm^3/atom
  [../]
  [./gamma]
    type = ParsedMaterial
    material_property_names = 'length_scale eVpJ T '
    # constant_names = 'gamma1 gamma0'
    # constant_expressions = '-4.5292e-04 2.5359'
    # function = '(gamma1*T + gamma0)*eVpJ*length_scale^2'
    constant_names = 'gamma_Tm'
    constant_expressions = '0.406'  # J/m^2
    function = 'gamma_Tm*eVpJ*length_scale^2'
    f_name = gamma
    outputs = exodus
  [../]

  ####### Calculation of kinetic parameters
  [./diffusivity]
    type = ParsedMaterial
    material_property_names = 'length_scale time_scale Em kb T D0'
    function = 'D0/length_scale^2*time_scale*exp(-Em/(kb*T))'
    f_name = D
    outputs = exodus
  [../]
  [./Mobility]
    type = ParsedMaterial
    f_name = M
    # material_property_names = 'kb T Va D'
    # args = cv
    # function = '1e=5*D/(kb*T)/Va'  ### Why 1e-5?
    material_property_names = 'D h fVcc:=D[fv,cvV,cvV] fBcc:=D[fb,cvB,cvB]'
    function = 'fcc:=fVcc*fBcc/((1-h)*fVcc+h*fBcc + 1e-30);D/fcc'
    outputs = exodus
  [../]

  ####### Calculation of free energy expressions
  [./fb_fit]
    type = DerivativeParsedMaterial
    args = 'cvB'
    material_property_names = 'Ab Bb eps'
    function = 'Ab*(sqrt(cvB^2 + eps^2) - eps) + 0.5*Bb*cvB^2'
    # f_name = fb_chem
    f_name = fb
  [../]
  [./fv_poly]
    type = DerivativeParsedMaterial
    args = cvV
    material_property_names = 'Av cVeq'
    function = '1/2 * Av * (cvV - cVeq)^2'
    f_name = fv
  [../]
  ####### Other functions needed for phase field method
  [./h]
    type = SwitchingFunctionMaterial
    function_name = h
    h_order = HIGH
    eta = psi
  [../]
  [./g]
    type = BarrierFunctionMaterial
    function_name = g
    g_order = SIMPLE
    eta = psi
  [../]

  ####### Define Phase Field Parameters
  [./kappa_eq]
    type = ParsedMaterial
    material_property_names = 'gamma delta'
    function = '3*sqrt(2)*gamma*delta'
    outputs = exodus
    f_name = kappa
  [../]

  ####### Values needed for quantitative nucleation
  [./nucleation_psi]
    type = DiscreteNucleation
    f_name = fnpsi
    op_names  = psi
    op_values = 1.0
    penalty = 88.2714
    penalty_mode = MIN
    map = mappsi
    outputs = exodus
  [../]

  [./prob_homeg]
    type = ParsedMaterial
    f_name = P_homog
    args = 'cv psi cv_avg'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    # function = 'if(psi>0.01,0,1)*if(cv<0,0,1)*cv*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    function = 'if(psi>0.01,0,1)*if(cv_avg<0,0,1)*cv_avg*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    # function = '0.0001'
    outputs = exodus
  [../]

  [./prob_heterog]
    type = ParsedMaterial
    f_name = P_heterog
    args = 'cv psi cv_avg'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em S'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    # function = 'if(psi>0.01,0,1)*if(cv<0,0,1)*cv*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
    function = 'if(psi>0.01,0,1)*if(cv_avg<0,0,1)*cv_avg*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
    # function = '0.0001'
    outputs = exodus
  [../]

  [./probability]
    type = ParsedMaterial
    f_name = P
    args = 'cv psi surf'
    material_property_names = 'P_homog P_heterog'
    function = '(surf*P_heterog + (1-surf)*P_homog)'
    outputs = exodus
  [../]

  [./shape_factor]
    type = ParsedMaterial
    f_name = S
    constant_names =       'pi       theta'
    constant_expressions = '3.141592 60'
    function = '(2+cos(theta*pi/180))*(1-(cos(theta*pi/180))^2)/4'
  [../]

  [./nucl_site_fraction]
    type = ParsedMaterial
    f_name = Ns
    material_property_names = 'Va'
    function = '1/Va'
    outputs = exodus
  [../]

  [./dgstar]
    type = ParsedMaterial
    material_property_names = 'pi gamma dgv'
    function = 'pi*gamma^2/(abs(dgv)+1e-30)'  #in 2D
    # function = '16*pi*gamma^3/3/(abs(dgv))^2'  #in 3D
    f_name = dgstar
    outputs = exodus
  [../]

  [./rstar]
    type = ParsedMaterial
    material_property_names = 'gamma dgv'
    function = 'gamma/(abs(dgv)+1e-30)'   # in 2D
    # function = '2*gamma/abs(dgv)'   # in 3D
    outputs = exodus
    f_name = rstar
  [../]

  [./dgvmatrix]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb eps'
    function = 'Ab*(sqrt((if(cv_avg<0,0,cv_avg)-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(if(cv_avg<0,0,cv_avg)-cBeq)^2'
    # f_name = fb_chem
    args = 'cv_avg'
    f_name = dgvmatrix
    outputs = exodus
  [../]

  [./dgvnucl]
    type = ParsedMaterial
    material_property_names = 'Av cVeq'
    args = 'cv_avg'
    # function = '1/2 * Av * (if(cv_avg<0,0,cv_avg) - cVeq)^2'
    function = '1/2 * Av * (cVeq - cVeq)^2' #???? This gives zero and you don't use cv_avg ????
    f_name = dgvnucl
    outputs = exodus
  [../]

  [./fcbulk_slope]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb eps'
    args = 'cv_avg'
    function = 'Ab*(if(cv_avg<0,0,cv_avg)/sqrt((if(cv_avg<0,0,cv_avg)-cBeq)^2+eps^2)) + Bb*if(cv_avg<0,0,cv_avg)'
    f_name = fcbulk_slope
    outputs = exodus
  [../]

  [./fcbulk_tangent]
    type = ParsedMaterial
    material_property_names = 'cVeq dgvmatrix fcbulk_slope'
    args = 'cv_avg'
    function = 'fcbulk_slope*(cVeq-if(cv_avg<0,0,cv_avg)) + dgvmatrix'
    f_name = fcbulk_tangent
  [../]

  [./dgv]
    type = ParsedMaterial
    material_property_names = 'dgvnucl fcbulk_tangent'
    function = 'dgvnucl - fcbulk_tangent'
    f_name = dgv
    outputs = exodus
  [../]

  [./cBeq_temp_dep]
    type = ParsedMaterial
    material_property_names = 'kb T Ef'
    constant_names =       'Sf'
    constant_expressions = '4.7' #???? What is Sf?????
    function = 'exp(Sf)*exp(-Ef/kb/T)'
    f_name = cBeq
    outputs = exodus
  [../]

  [./nucl_radius]
    type = ParsedMaterial
    function = '1'
    f_name = r_crit
  [../]


[]

[Kernels]
  # # enforce c = (1-h(eta))*cm + h(eta)*cd
  # [./PhaseConc]
  #   type = KKSPhaseConcentration
  #   variable = cvV
  #   c        = cv
  #   eta      = psi
  # [../]
  # # enforce pointwise equality of chemical potentials
  # [./ChemPotVacancies]
  #   type = KKSPhaseChemicalPotential
  #   variable = cvB
  # [../]

  # Cahn-Hilliard Equation
  #
  [./CHBulk]
    type = KKSCHBulk
    variable = cv
    mob_name = M
    args = 'psi'
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = cv
  [../]

  #
  # Allen-Cahn Equation
  #
  # [./ACBulkF]
  #   type = KKSACBulkF
  #   variable = psi
  #   args     = 'cv cvB cvV'
  #   w        = 88.27 #kappa/delta^2
  # [../]
  # [./ACBulkC]
  #   type = KKSACBulkC
  #   variable = psi
  #   args = 'cv cvB cvV'
  # [../]
  # [./ACInterface]
  #   type = ACInterface
  #   variable = psi
  #   kappa_name = kappa
  # [../]
  # [./detadt]
  #   type = TimeDerivative
  #   variable = psi
  # [../]
  #
  # [./fnuclpsi]
  #   type = AllenCahn
  #   variable = psi
  #   mob_name = M
  #   f_name = fnpsi
  # [../]

  # [./cv_nucl_force_test]
  #   type = DiscreteNucleationForce
  #   variable = cv
  #   map = mappsi
  # [../]

[]

[AuxKernels]
  [./f_dens]
    variable = f_dens
    type = KKSGlobalFreeEnergy
    w = 88.27 #kappa/delta^2
  [../]

  [./psi_nucl_map]
    type = DiscreteNucleationAux
    map = mappsi
    variable = psi_map
    no_nucleus_value = 0
    nucleus_value = 1
  [../]

  [./av_cv]
    type = FunctionAux
    variable = cv_avg
    function = get_av
  [../]
[]

[Functions]
  [./get_av]
    type = ParsedFunction
    vars = 'avg_c'
    vals = 'avg_c'
    value = avg_c
  [../]
[]

[BCs]
  [./left_flux]
    type = NeumannBC
    boundary = left
    variable = cv
    value = 1.2e-7 #Flux
  [../]
  # [./Periodic]
  #   [./All]
  #     auto_direction = 'x y'
  #   [../]
  # [../]
[]

[UserObjects]
  [./inserter]
    type = DiscreteNucleationInserter
    hold_time = 1e-2 #lower limit = 7.2e-5
    probability = P
    radius = r_crit
    # probability = 0.0001 # just for test
    time_dependent_statistics = False
  [../]

  [./mappsi]
    type = DiscreteNucleationMap
    # periodic = psi
    periodic = cv
    inserter = inserter
    int_width = 1
    # radius = nrad
  [../]

  [./term]
    type = Terminator
    expression = 'rate > 100'
  [../]
[]


[Postprocessors]
  [./total_F]
    type = ElementIntegralVariablePostprocessor
    variable = f_dens
  [../]
  [./total_cv]
    type = ElementIntegralVariablePostprocessor
    variable = cv
    execute_on = 'initial timestep_end'
  [../]
  [./avg_c]
    type = ElementAverageValue
    variable = cv
    execute_on = 'initial timestep_end'
  [../]
  [./porosity]
    type = ElementAverageValue
    variable = psi
    execute_on = 'initial timestep_end'
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./max_cv]
    type = ElementExtremeValue
    variable = cv
    value_type = max
    execute_on = 'initial timestep_end'
  [../]
  [./min_cv]
    type = ElementExtremeValue
    variable = cv
    value_type = min
    execute_on = 'initial timestep_end'
  [../]

  [./rate]
    type = DiscreteNucleationData
    value = RATE
    inserter = inserter
  [../]
  [./dtnuc]
    type = DiscreteNucleationTimeStep
    inserter = inserter
    p2nucleus = 0.0001
    dt_max = 10
  [../]
  [./update]
    type = DiscreteNucleationData
    value = UPDATE
    inserter = inserter
  [../]
  [./count]
    type = DiscreteNucleationData
    value = COUNT
    inserter = inserter
  [../]

[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  # solve_type = NEWTON
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'lu'
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'bjacobi  gmres'

  scheme = bdf2
  # end_time = 30

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  l_tol = 1e-05
  nl_max_its = 15
  l_max_its = 30
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    cutback_factor = 0.8
    growth_factor = 1.2
    # optimal_iterations = 7
    # linear_iteration_ratio = 100
    # timestep_limiting_postprocessor = dtnuc
  [../]

  # [./Adaptivity]
  #   max_h_level = 1
  #   initial_adaptivity = 1
  #   refine_fraction = 0.7
  #   coarsen_fraction = 0.05
  # [../]
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]

[Debug]
  show_var_residual_norms = true
[]

# [MultiApps]
#   [./ECPsolver]
#     type = FullSolveMultiApp
#     app_type = discrete_nuclApp
#     execute_on = 'FINAL'
#     positions = '0 0 0'
#     input_files = ./ecp_solver_initial.i
#   [../]
# []
#
# [Transfers]
#   [./Get_cvV_ECPsolver_to_PF]
#     type = MultiAppMeshFunctionTransfer
#     direction = from_multiapp
#     multi_app = ECPsolver
#     source_variable = cvV
#     variable = cvV
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Get_cvB_ECPsolver_to_PF]
#     type = MultiAppMeshFunctionTransfer
#     direction = from_multiapp
#     multi_app = ECPsolver
#     source_variable = cvB
#     variable = cvB
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_psi_PF_to_ECPsolver]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = ECPsolver
#     source_variable = psi
#     variable = psi
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_cv_PF_to_ECPsolver]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = ECPsolver
#     source_variable = cv
#     variable = cv
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_cvV_PF_to_ECPsolver]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = ECPsolver
#     source_variable = cvV
#     variable = cvV
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_cvB_PF_to_ECPsolver]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = ECPsolver
#     source_variable = cv_avg
#     variable = cvB
#     execute_on = SAME_AS_MULTIAPP
#   [../]
# []
