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
  [./surf]
    order = FIRST
    family = LAGRANGE
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
    value = 1e-4
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
    value = 1e-4
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

  # [./cv_circle_IC]
  #   type = SmoothCircleIC
  #   variable = cv
  #   radius = 2
  #   int_width = 1
  #   profile = TANH
  #   invalue = 0.999
  #   outvalue = 0.01
  #   x1 = 50
  #   y1 = 50
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
  [./surface]
    type = ParsedFunction
    vars = 'xl yl'
    vals = '2  100'
    value = 'if(x<xl,1,0)'
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
    prop_values = '3.6 300 1e-9         1e-9       6.24150934e+18 8.6173303e-5  0.6928  0.055 0.52 270.5 32 200 1'
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
    args = 'cv psi'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(psi>0.01,0,1)*cv*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    # function = '0.0001'
    outputs = exodus
  [../]

  [./prob_heterog]
    type = ParsedMaterial
    f_name = P_heterog
    args = 'cv psi'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em S'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(psi>0.01,0,1)*cv*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
    # function = '0.0001'
    outputs = exodus
  [../]

  [./probability]
    type = ParsedMaterial
    f_name = P
    args = 'cv psi surf'
    material_property_names = 'P_homog P_heterog'
    function = 'surf*P_heterog + (1-surf)*P_homog'
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

  [./tau]
    type = ParsedMaterial
    f_name = tau
    material_property_names = 'D0 cBeq cVeq'
    constant_names =       'cinit radius'
    constant_expressions = '1e-4 0.14'
    function = 'radius^2/D0*(cVeq-cinit)^2/(cinit-cBeq)^2'
    outputs = exodus
  [../]

  [./dgstar]
    type = ParsedMaterial
    constant_names = 'pi'
    constant_expressions = '3.141592'
    material_property_names = 'gamma dgv'
    function = 'pi*gamma^2/abs(dgv)'  #in 2D
    # function = '16*pi*gamma^3/3/(abs(dgv))^2'  #in 3D
    f_name = dgstar
    outputs = exodus
  [../]

  [./rstar]
    type = ParsedMaterial
    material_property_names = 'gamma dgv'
    function = 'gamma/abs(dgv)'   # in 2D
    # function = '2*gamma/abs(dgv)'   # in 3D
    outputs = exodus
    f_name = rstar
  [../]

  [./dgvmatrix]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb'
    constant_names =       'cinit eps'
    constant_expressions = '1e-4 0.001'
    function = 'Ab*(sqrt((cinit-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(cinit-cBeq)^2'
    # f_name = fb_chem
    f_name = dgvmatrix
    outputs = exodus
  [../]
  [./dgvnucl]
    type = ParsedMaterial
    material_property_names = 'Av cVeq'
    constant_names =       'cinput'
    constant_expressions = '1.0'
    function = '1/2 * Av * (cinput - cVeq)^2'
    f_name = dgvnucl
    outputs = exodus
  [../]

  [./fcbulk_slope]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb'
    constant_names =       'cinit eps'
    constant_expressions = '1e-4 0.001'
    function = 'Ab*(cinit/sqrt((cinit-cBeq)^2+eps^2)) + Bb*cinit'
    f_name = fcbulk_slope
    outputs = exodus
  [../]

  [./fcbulk_tangent]
    type = ParsedMaterial
    material_property_names = 'cVeq dgvmatrix fcbulk_slope'
    constant_names =       'cinit'
    constant_expressions = '1e-4'
    function = 'fcbulk_slope*(cVeq-cinit) + dgvmatrix'
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
    constant_expressions = '4.7'
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

  [./fnuclpsi]
    type = AllenCahn
    variable = psi
    mob_name = M
    f_name = fnpsi
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

  [./construct_surface]
    type = FunctionAux
    function = surface
    variable = surf
  [../]
[]

[BCs]
  [./left_flux]
    type = NeumannBC
    boundary = left
    variable = cv
    value = 1e-2 #Flux
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
    hold_time = 1e-4 #lower limit = 7.2e-5
    probability = P
    radius = r_crit
    # probability = 0.0001 # just for test
    # time_dependent_statistics = False
  [../]

  [./mappsi]
    type = DiscreteNucleationMap
    periodic = psi
    inserter = inserter
    int_width = 1
    # radius = nrad
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
  end_time = 30

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  l_tol = 1e-05
  nl_max_its = 15
  l_max_its = 30
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    cutback_factor = 0.8
    growth_factor = 1.2
    # optimal_iterations = 7
    # linear_iteration_ratio = 100
    timestep_limiting_postprocessor = dtnuc
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

[Debug]
  show_var_residual_norms = true
[]
