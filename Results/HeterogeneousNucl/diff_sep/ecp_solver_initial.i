#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
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

  # radius = 0.10
  # int_width = 0.04
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

  # order parameter
  [./psi]
  [../]
  # vacancy concentration
  [./cv]
  [../]
[]

[Variables]

  # bulk phase concentration (matrix)
  [./cvB]
  [../]
  # hydrogen phase concentration (delta phase)
  [./cvV]
  [../]
[]

[ICs]
  # [./cv_IC]
  #   type = ConstantIC
  #   variable = cv
  #   value = 2.021e-7
  # [../]
  # [./psi_IC]
  #   type = ConstantIC
  #   variable = psi
  #   value = 0
  # [../]
  # [./cvV_IC]
  #   type = ConstantIC
  #   variable = cvV
  #   value = 1.0
  # [../]
  # [./cvB_IC]
  #   type = ConstantIC
  #   variable = cvB
  #   value = 2.021e-7
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
    constant_expressions = '0.01 270.5 32'
    function = 'A*(sqrt(cvB^2 + eps^2) - eps) + 0.5*B*cvB^2'
    # f_name = fb_chem
    f_name = fb
  [../]
  [./fv_poly]
    type = DerivativeParsedMaterial
    args = cvV
    constant_names = 'A ceq'
    constant_expressions = '1000 1.0'
    function = '1/2 * A * (cvV - ceq)^2'
    f_name = fv

    # type = DerivativeParsedMaterial
    # args = 'cvV'
    # constant_names =       'eps   A     B'
    # constant_expressions = '0.01 270.5 32'
    # function = 'A*(sqrt((cvV-1)^2 + eps^2) - eps) + 0.5*B*(cvV-1)^2'
    # # f_name = fb_chem
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
    constant_names = 'delta'
    constant_expressions = '0.6928'
    material_property_names = 'gamma'
    function = '3*sqrt(2)*gamma*delta'
    outputs = exodus
    f_name = kappa
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
  [./gamma_surf]
    type = ParsedMaterial
    material_property_names = 'length_scale eVpJ T '
    # constant_names = 'gamma1 gamma0'
    # constant_expressions = '-4.5292e-04 2.5359'
    # function = '(gamma1*T + gamma0)*eVpJ*length_scale^2'
    constant_names = 'gamma'
    constant_expressions = '0.406'  # J/m^2
    function = 'gamma*eVpJ*length_scale^2'
    f_name = gamma_surf
    outputs = exodus
  [../]
  [./kappa_L]
    type = GenericConstantMaterial
    # prop_names = 'kappa L'
    # prop_values = '24 10.0'
    prop_names =  'L   T   length_scale time_scale eVpJ           kb            delta   Em    Ef   Ab    Bb Av  cVeq'
    prop_values = '3.6 300 1e-9         1e-9       6.24150934e+18 8.6173303e-5  0.6928  0.055 0.52 270.5 32 1000 1'
  [../]
  [./Va]
    type = ParsedMaterial
    f_name = Va
    # material_property_names = 'length_scale'
    # function = '5.586e24*length_scale^3'
    function = '2.024e-4'  #nm^3/atom
  [../]

  [./Mobility]
    type = ParsedMaterial
    f_name = M
    material_property_names = 'kb T Va D'
    args = cv
    function = '1e-5*D/(kb*T)/Va'
    outputs = exodus
  [../]

  [./diffusivity]
    type = ParsedMaterial
    material_property_names = 'length_scale time_scale Em kb T'
    constant_names = 'D0'
    constant_expressions = '8.37707e-7'  # m^2/s
    # function = '1.39*(1e-2/length_scale^2)*time_scale*exp(-Em/(kb*T))'
    function = 'D0/length_scale^2*time_scale*exp(-Em/(kb*T))'
    f_name = D
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
  #
  # Cahn-Hilliard Equation
  #
  # [./CHBulk]
  #   type = KKSCHBulk
  #   variable = cv
  #   mob_name = M
  #   args = 'psi'
  # # [../]
  # [./dcdt]
  #   type = TimeDerivative
  #   variable = cv
  # [../]

  #
  # Allen-Cahn Equation
  #
  # [./ACBulkF]
  #   type = KKSACBulkF
  #   variable = psi
  #   args     = 'cv cvB cvV'
  #   # w        = 88.27 #kappa/delta^2
  #   w = 2.238380e+01
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

[]

[AuxKernels]
  [./f_dens]
    variable = f_dens
    type = KKSGlobalFreeEnergy
    # w = 88.27 #kappa/delta^2
    w = 88.2714
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

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  # solve_type = NEWTON
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'asm'
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'bjacobi  gmres'

  # scheme = bdf2
  # end_time = 30
  #
  # nl_rel_tol = 1e-9
  # nl_abs_tol = 1e-9
  # l_tol = 1e-05
  nl_max_its = 100
  l_max_its = 100
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

[Debug]
  show_var_residual_norms = true
[]
