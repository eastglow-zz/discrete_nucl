#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 100
  ymax = 100
[]

[GlobalParams]
  fa_name  = fb
  fb_name  = fv
  ca       = cvB
  cb       = cvV
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
[]

[Variables]
  # order parameter
  [./psi]
  [../]
  # vacancy concentration
  [./cv]
  [../]
  # bulk phase concentration (matrix)
  [./cvB]
  [../]
  # hydrogen phase concentration (delta phase)
  [./cvV]
  [../]
[]

[ICs]
  # [./cv_circle_IC]
  #  type = LatticeSmoothCircleIC
  #  variable = cv
  #  circles_per_side = '4 4'
  #  radius = 2
  #  pos_variation = 8.0
  #  int_width = 1
  #  profile = TANH
  #  invalue = 0.999
  #  outvalue = 1e-9
  # [../]
  #
  # [./psi_circle_IC]
  #  type = LatticeSmoothCircleIC
  #  variable = psi
  #  circles_per_side = '4 4'
  #  radius = 2
  #  pos_variation = 8.0
  #  int_width = 2
  #  profile = TANH
  #  invalue = 1
  #  outvalue = 0
  # [../]

  [./cv_IC]
    type = ConstantIC
    variable = cv
    value = 1e-9
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
    value = 1e-9
  [../]
[]

[Materials]
  # [./fb]
  #   type = DerivativeParsedMaterial
  #   args = 'cvB'
  #   constant_names = 'length_scale Ef kb T Va ev tol'
  #   constant_expressions = '1e-9 1.6 8.6173303e-5 500 5.586e24*length_scale^3 Ef/Va 1e-8' #For alpha-Fe
  #   function = 'cvB * ev + kb * T / Va * (cvB * plog(cvB,tol) + (1 - cvB) * plog(1 - cvB,tol))'
  #   # f_name = fb_chem
  #   f_name = fb
  # [../]

  [./fb_fit]
    type = DerivativeParsedMaterial
    args = 'cvB'
    constant_names =       'eps   A   B'
    constant_expressions = '0.001 260 50'
    function = 'A*(sqrt(cvB^2 + eps^2) - eps) + 0.5*B*cvB^2'
    # f_name = fb_chem
    f_name = fb
  [../]
  [./fv_poly]
    type = DerivativeParsedMaterial
    args = cvV
    constant_names = 'A ceq'
    constant_expressions = '200 1.0'
    function = '1/2 * A * (cvV - ceq)^2'
    f_name = fv
  [../]

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
  [./kappa_eq]
    type = ParsedMaterial
    constant_names = 'length_scale eVpJ T gamma1 gamma0 delta'
    constant_expressions = '1e-9 6.24150934e+18 500 -4.5292e-04 2.5359 0.6928'
    function = 'gamma:=(gamma1*T + gamma0)*eVpJ*length_scale^2; 3*sqrt(2)*gamma*delta'
    outputs = exodus
    f_name = kappa
  [../]
  [./kappa_L]
    type = GenericConstantMaterial
    # prop_names = 'kappa L'
    # prop_values = '24 10.0'
    prop_names =  'L   One'
    prop_values = '3.6 1.0'
  [../]
  [./Mobility]
    type = ParsedMaterial
    f_name = M
    constant_names = 'length_scale time_scale Em kb T D0 Va'
    constant_expressions = '1e-9 1e-6 0.89 8.6173303e-5 500 1.39*(1e-2/length_scale^2)*time_scale 5.586e24*length_scale^3'
    args = cv
    function = '1e-5*D0*exp(-Em/(kb*T))/(kb*T)/Va'
    outputs = exodus
  [../]

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

  [./probability]
    type = ParsedMaterial
    f_name = P
    args = cv
    function = cv
    outputs = exodus
  [../]

[]

[Kernels]
  # enforce c = (1-h(eta))*cm + h(eta)*cd
  [./PhaseConc]
    type = KKSPhaseConcentration
    variable = cvV
    c        = cv
    eta      = psi
  [../]
  # enforce pointwise equality of chemical potentials
  [./ChemPotVacancies]
    type = KKSPhaseChemicalPotential
    variable = cvB
  [../]
  #
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
  [./ACBulkF]
    type = KKSACBulkF
    variable = psi
    args     = 'cv cvB cvV'
    w        = 88.27 #kappa/delta^2
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = psi
    args = 'cv cvB cvV'
  [../]
  [./ACInterface]
    type = ACInterface
    variable = psi
    kappa_name = kappa
  [../]
  [./detadt]
    type = TimeDerivative
    variable = psi
  [../]

  # [./fnuclpsi]
  #   type = AllenCahn
  #   variable = psi
  #   mob_name = M
  #   f_name = fnpsi
  # [../]

  [./psi_nucl_force]
    type = DiscreteNucleationForce
    variable = psi
    map = mappsi
    no_nucleus_value = 0
    nucleus_value = 1
  [../]
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
[]

[BCs]
  [./left_flux]
    type = NeumannBC
    boundary = left
    variable = cv
    value = 0.8e-1 #Flux
  [../]
[]

[UserObjects]
  [./inserter]
    type = DiscreteNucleationInserter
    hold_time = 1e-5
    # probability = P
    probability = 1 # just for test
  [../]

  [./mappsi]
    type = DiscreteNucleationMap
    radius = 2
    periodic = psi
    inserter = inserter
    int_width = 1
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
  end_time = 30

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  l_tol = 1e-04
  nl_max_its = 30
  l_max_its = 100
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.85e-4
    cutback_factor = 0.75
    growth_factor = 1.1
    optimal_iterations = 7
    linear_iteration_ratio = 100
  [../]
  [./Adaptivity]
    max_h_level = 1
    initial_adaptivity = 1
    refine_fraction = 0.7
    coarsen_fraction = 0.05
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
