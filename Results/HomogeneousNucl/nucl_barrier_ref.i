#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  # xmax = 100
  # ymax = 100
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

  [./f_inter]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./f_bulk]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./grad_x_phi]
    order = FIRST
    family = MONOMIAL
  [../]

  [./grad_y_phi]
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

  [./cv_IC]
    type = ConstantIC
    variable = cv
    value = 0.001
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
    value = 0.001
  [../]

  # [./psi_circle_IC]
  #   type = SmoothCircleIC
  #   variable = psi
  #   radius = 2
  #   int_width = 0.6928
  #   profile = TANH
  #   invalue = 1
  #   outvalue = 0
  #   x1 = 25
  #   y1 = 25
  # [../]
  #
  # [./cv_circle_IC]
  #   type = SmoothCircleIC
  #   variable = cv
  #   radius = 2
  #   int_width = 0.6928
  #   profile = TANH
  #   invalue = 1
  #   outvalue = 0.001
  #   x1 = 25
  #   y1 = 25
  # [../]

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
    vals =  '50 50 0    1    0.001  1  31.6386'
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
    constant_expressions = '0.001 270.5 32'
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
    constant_names = 'delta'
    constant_expressions = '0.6928'
    material_property_names = 'gamma'
    function = '3*sqrt(2)*gamma*delta'
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

  [./diffusivity]
    type = ParsedMaterial
    constant_names = 'length_scale time_scale'
    constant_expressions = '1e-9 1e-6 '
    function = '1.39*(1e-2/length_scale^2)*time_scale'
    f_name = Dop
  [../]

  [./Inter_f_density]
    type = ParsedMaterial
    material_property_names = 'kappa'
    args = 'grad_x_phi grad_y_phi'
    function = 'kappa*(grad_x_phi^2 + grad_y_phi^2)'
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

[]

[AuxKernels]
  [./f_dens]
    variable = f_dens
    type = KKSGlobalFreeEnergy
    w = 88.27 #kappa/delta^2
  [../]

  [./f_dens_inter]
    variable = f_inter
    type = MaterialRealAux
    property = f_dens_inter
    execute_on = LINEAR
  [../]

  [./f_dens_bulk]
    variable = f_bulk
    type = MaterialRealAux
    property = f_dens_bulk
    execute_on = LINEAR
  [../]

  [./dphidx]
    type = VariableGradientComponent
    gradient_variable = grad_x_phi
    variable = phi
    component = x
  [../]

  [./dphidy]
    type = VariableGradientComponent
    gradient_variable = grad_y_phi
    variable = phi
    component = y
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
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
  l_tol = 1e-05
  nl_max_its = 15
  l_max_its = 30
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.85e-3
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
