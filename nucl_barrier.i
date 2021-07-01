#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 200
  ny = 200
  # xmax = 100
  # ymax = 100
  xmax = 10
  ymax = 10
[]

[GlobalParams]
  fa_name  = fb
  fb_name  = fv
  ca       = cvB
  cb       = cvV

  radius = 0.2
  int_width = 0.1

[]

[AuxVariables]
  [./f_dens]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./f_inter]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./f_bulk]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./grad_x_psi]
    order = FIRST
    family = MONOMIAL
  [../]

  [./grad_y_psi]
    order = FIRST
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

  # [./cv_IC]
  #   type = ConstantIC
  #   variable = cv
  #   value = 0.001
  # [../]
  # [./psi_IC]
  #   type = ConstantIC
  #   variable = psi
  #   value = 0
  # [../]
  [./cvV_IC]
    type = SmoothCircleIC
    variable = cvV
    radius = 0.2
    int_width = 0.1
    profile = TANH
    invalue = 1.0
    outvalue = 1.0
    x1 = 0
    y1 = 0
  [../]
  [./cvB_IC]
    type = SmoothCircleIC
    profile = TANH
    variable = cvB
    invalue = 0.0001
    outvalue = 0.0001
    x1 = 0
    y1 = 0
  [../]

  [./psi_circle_IC]
    type = SmoothCircleIC
    variable = psi
    profile = TANH
    invalue = 1
    outvalue = 0
    x1 = 0
    y1 = 0
  [../]

  [./cv_circle_IC]
    type = SmoothCircleIC
    variable = cv
    profile = TANH
    invalue = 1.0
    outvalue = 0.0001
    x1 = 0
    y1 = 0
  [../]

  # [./cv_circle_IC]
  #   type = FunctionIC
  #   variable = cv
  #   function = cv_profile
  # [../]
[]

[Functions]
  [./cv_profile]
    type = ParsedFunction
    vars =  'x0 y0 cBeq cNeq cin    rN rd'
    vals =  '25 25 0    1    0.001  1  31.6386'
    value = 'r:=sqrt((x-x0)^2+(y-y0)^2);if(r <= rN, cNeq, if(r <= rN+rd, cBeq + (cin - cBeq)*r/rd, cin))'
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
    constant_names =       'eps   A     B'
    constant_expressions = '0.01 270.5 32'
    function = 'A*(sqrt(cvB^2 + eps^2) - eps) + 0.5*B*cvB^2'
    # f_name = fb_chem
    f_name = fb
  [../]
  [./fv_poly]
    type = DerivativeParsedMaterial
    args = cvV
    constant_names =       'A   ceq'
    constant_expressions = '1000 1.0'
    function = '1/2 * A * (cvV - ceq)^2'
    f_name = fv

    # type = DerivativeParsedMaterial
    # args = 'cvV'
    # constant_names =       'eps   A     B'
    # constant_expressions = '0.01 270.5 32'
    # function = 'A*(sqrt((cvV-1)^2 + eps^2) - eps) + 0.5*B*(cvV-1)^2'
    # f_name = fv
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
  # [./kappa_eq]
  #   type = ParsedMaterial
  #   constant_names = 'length_scale eVpJ T gamma1 gamma0 delta'
  #   constant_expressions = '1e-9 6.24150934e+18 500 -4.5292e-04 2.5359 0.6928'
  #   function = 'gamma:=(gamma1*T + gamma0)*eVpJ*length_scale^2; 3*sqrt(2)*gamma*delta'
  #   outputs = exodus
  #   f_name = kappa
  # [../]
  [./kappa_eq]
    type = ParsedMaterial
    material_property_names = 'gamma delta alpha'
    # function = '3*sqrt(2)*gamma*delta'
    function = '3*gamma*delta/alpha'
    outputs = exodus
    f_name = kappa
  [../]
  [./gamma]
    type = ParsedMaterial
    constant_names = 'length_scale eVpJ T gamma1 gamma0'
    constant_expressions = '1e-9 6.24150934e+18 500 -4.5292e-04 2.5359'
    # function = '(gamma1*T + gamma0)*eVpJ*length_scale^2'
    function = '0.406*eVpJ*length_scale^2'
    f_name = gamma
    outputs = exodus
  [../]
  [./consts]
    type = GenericConstantMaterial
    # prop_names = 'kappa L'
    # prop_values = '24 10.0'
    prop_names =  'One  delta   alpha'
    prop_values = '1.0  0.1       2.9444' # alpha = ln(0.95/0.05), interface is defined where psi ranged [0.05, 0.95]
  [../]
  [./Mobility]
    type = ParsedMaterial
    f_name = M
    constant_names = 'length_scale time_scale Em kb T D0 Va'
    constant_expressions = '1e-9 1e-6 0.89 8.6173303e-5 500 1.39*(1e-2/length_scale^2)*time_scale 5.586e24*length_scale^3'
    args = cv
    function = '1e-5*D0*exp(-Em/(kb*T))/(kb*T)/Va'  #original
    # function = '1e-4*D0*exp(-Em/(kb*T))/(kb*T)/Va'
    outputs = exodus
  [../]

  [./diffusivity]
    type = ParsedMaterial
    constant_names = 'length_scale time_scale'
    constant_expressions = '1e-9 1e-6 '
    function = '1.39*(1e-2/length_scale^2)*time_scale'
    f_name = Dop
  [../]

  [./pf_mobility]
    type = ParsedMaterial
    f_name = L
    material_property_names = 'M gamma delta'
    function = '0.01'
  [../]

  [./Inter_f_density]
    type = ParsedMaterial
    material_property_names = 'kappa'
    args = 'grad_x_psi grad_y_psi'
    function = 'kappa*(grad_x_psi^2 + grad_y_psi^2)'
    f_name = f_dens_inter
    outputs = exodus
    execute_on = 'LINEAR INITIAL'
  [../]

  [./fbulk_density]
    type = ParsedMaterial
    material_property_names = 'f_dens_inter'
    args = 'f_dens'
    f_name = f_dens_bulk
    function = 'f_dens - f_dens_inter'
    outputs = exodus
    execute_on = 'LINEAR INITIAL'
  [../]

  [./w_barrier]
    type = ParsedMaterial
    material_property_names = 'gamma delta alpha'
    f_name = doublewell
    function = '3*alpha*gamma/delta'
    execute_on = 'LINEAR INITIAL'
    outputs = exodus
  [../]

[]

[Kernels]
  # enforce c = (1-h(eta))*cm + h(eta)*cd
  [./PhaseConcV]
    type = KKSPhaseConcentration
    variable = cvV
    c        = cv
    eta      = psi
  [../]
  # enforce pointwise equality of chemical potentials
  [./ChemPotB]
    type = KKSPhaseChemicalPotential
    variable = cvB
  [../]

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
    # w        = 88.27 #kappa/delta^2
    w = 2.238380e+01
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

[]

[AuxKernels]
  [./f_dens]
    variable = f_dens
    type = KKSGlobalFreeEnergy
    # w = 88.27 #kappa/delta^2
    w = 2.238380e+01
    execute_on = 'LINEAR INITIAL'
  [../]

  [./f_dens_inter]
    variable = f_inter
    type = MaterialRealAux
    property = f_dens_inter
    execute_on = 'LINEAR INITIAL'
  [../]

  [./f_dens_bulk]
    variable = f_bulk
    type = MaterialRealAux
    property = f_dens_bulk
    execute_on = 'LINEAR INITIAL'
  [../]

  [./dphidx]
    type = VariableGradientComponent
    gradient_variable = psi
    variable = grad_x_psi
    component = x
    execute_on = 'LINEAR INITIAL'
  [../]

  [./dphidy]
    type = VariableGradientComponent
    gradient_variable = psi
    variable = grad_y_psi
    component = y
    execute_on = 'LINEAR INITIAL'
  [../]
[]

[BCs]
  # [./Periodic]
  #   [./All]
  #     auto_direction = 'x y'
  #   [../]
  # [../]
[]

[Postprocessors]
  [./total_F]
    type = ElementIntegralVariablePostprocessor
    variable = f_dens
  [../]
  [./inter_F]
    type = ElementIntegralVariablePostprocessor
    variable = f_inter
  [../]
  [./bulk_F]
    type = ElementIntegralVariablePostprocessor
    variable = f_bulk
  [../]

  [./kappa_out]
    type = ElementAverageMaterialProperty
    mat_prop = kappa
  [../]

  [./w_out]
    type = ElementAverageMaterialProperty
    mat_prop = doublewell
  [../]


  # [./dt]
  #   type = TimestepSize
  # [../]
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
  l_tol = 1e-05
  nl_max_its = 15
  l_max_its = 30
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
    cutback_factor = 0.8
    growth_factor = 1.2
    optimal_iterations = 7
    linear_iteration_ratio = 100
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


[MultiApps]
  [./ECPsolver]
    type = FullSolveMultiApp
    app_type = discrete_nuclApp
    execute_on = 'INITIAL'
    positions = '0 0 0'
    input_files = ./ecp_solver_initial.i
  [../]
[]

[Transfers]
  [./Get_cvV_ECPsolver_to_PF]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = ECPsolver
    source_variable = cvV
    variable = cvV
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Get_cvB_ECPsolver_to_PF]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = ECPsolver
    source_variable = cvB
    variable = cvB
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_psi_PF_to_ECPsolver]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = ECPsolver
    source_variable = psi
    variable = psi
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_cv_PF_to_ECPsolver]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = ECPsolver
    source_variable = cv
    variable = cv
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_cvV_PF_to_ECPsolver]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = ECPsolver
    source_variable = cvV
    variable = cvV
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_cvB_PF_to_ECPsolver]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = ECPsolver
    source_variable = cvB
    variable = cvB
    execute_on = SAME_AS_MULTIAPP
  [../]
[]
