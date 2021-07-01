#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50 #Number of elements in the x-direction
  ny = 50 #Number of elements in the y-direction
  xmax = 100 #X-direction domain size, for this file it is in nm
  ymax = 100 #Y-direction domain size, for this file it is in nm
  uniform_refine = 2
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

  [./cv_analytical]
  [../]

  [./time]
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

[Functions]
  [./get_av]
    type = ParsedFunction
    vars = 'avg_c'
    vals = 'avg_c'
    value = avg_c
  [../]
  [./Left_BC]
    type = ParsedFunction
    #vars = flux (mol/m^2-s) length_scale time_scale Va (nm^2/atom) Na (atoms/mol)
    #The length_scale, time_scale and Va here should match in materials block
    vars = 'flux length_scale time_scale Va Na         D'
    vals = '1e-2 length_scale time_scale Va 6.02214e23 diffusivity_avg'
    value = 'flux*Na*length_scale^2*time_scale*Va/D'
  [../]
  [./Left_slope]
    type = ParsedFunction
    #vars = flux (mol/m^2-s) length_scale time_scale Va (nm^3/atom) Na (atoms/mol)
    #The length_scale, time_scale and Va here should match in materials block
    vars = 'flux length_scale time_scale Va Na         D'
    vals = '1e-2 length_scale time_scale Va 6.02214e23 diffusivity_avg'
    value = 'flux*Na*length_scale^2*time_scale*Va/D'
  [../]

  [./analytical]
    type = ParsedFunction
    vars = 'D dCdx C0'
    vals = 'diffusivity_avg Left_slope 2e-7'
    value = 'C0 + dCdx*(2*sqrt(D*t/pi)*exp(-x^2/(4*D*t)) - x*(1-erf(x/(2*sqrt(D*t)))))'
  [../]

  [./no_gradient]
    type = ParsedFunction
    vars = 'D dCdx C0 lx'
    vals = 'diffusivity_avg Left_slope 2e-7 100'
    value = 'C0 + dCdx*D/lx*t'
  [../]

  [./time_function]
    type = ParsedFunction
    value = 't'
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

  [./Va] #Atomic volume
    type = ParsedMaterial
    f_name = Va
    material_property_names = 'length_scale'
    constant_names = 'Va'
    constant_expressions = 20.24 #In angstroms^3/atom (Ma and Dudarev 2019)
    function = 'Va*(1e-10/length_scale)^3'  #Change units of length
  [../]

  [./gamma]
    type = ParsedMaterial
    material_property_names = 'length_scale eVpJ T '
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
    args = 'cv psi cv_analytical'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(psi>0.01,0,1)*if(cv_analytical<0,0,1)*cv_analytical*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    outputs = exodus
  [../]

  [./prob_heterog]
    type = ParsedMaterial
    f_name = P_heterog
    args = 'cv psi cv_analytical'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em S'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(psi>0.01,0,1)*if(cv_analytical<0,0,1)*cv_analytical*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
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
    function = 'Ab*(sqrt((if(cv_analytical<0,0,cv_analytical)-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(if(cv_analytical<0,0,cv_analytical)-cBeq)^2'
    # f_name = fb_chem
    args = 'cv_analytical'
    f_name = dgvmatrix
    outputs = exodus
  [../]

  [./dgvnucl]
    type = ParsedMaterial
    material_property_names = 'Av cVeq'
    args = 'cv_analytical'
    function = '1/2 * Av * (cVeq - cVeq)^2' #This sets the energy at nucleation, and currently is set to zero
    f_name = dgvnucl
    outputs = exodus
  [../]

  [./fcbulk_slope]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb eps'
    args = 'cv_analytical'
    function = 'Ab*(if(cv_analytical<0,0,cv_analytical)/sqrt((if(cv_analytical<0,0,cv_analytical)-cBeq)^2+eps^2)) + Bb*if(cv_analytical<0,0,cv_analytical)'
    f_name = fcbulk_slope
    outputs = exodus
  [../]

  [./fcbulk_tangent]
    type = ParsedMaterial
    material_property_names = 'cVeq dgvmatrix fcbulk_slope'
    args = 'cv_analytical'
    function = 'fcbulk_slope*(cVeq-if(cv_analytical<0,0,cv_analytical)) + dgvmatrix'
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
  # Cahn-Hilliard Equation
  #

  # [./diffusion]
  #   type = MatDiffusion
  #   variable = cv
  #   diffusivity = D
  # [../]
  [./dcdt]
    type = TimeDerivative
    variable = cv
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

  [./av_cv]
    type = FunctionAux
    variable = cv_avg
    function = get_av
  [../]

  [./cv_analytic]
    type = FunctionAux
    variable = cv_analytical
    function = analytical
    # function = no_gradient
  [../]

  [./time]
    type = FunctionAux
    variable = time
    function = time_function
  [../]
[]

[BCs]
  [./left_flux]
    type = FunctionNeumannBC
    boundary = left
    variable = cv
    function = Left_BC
  [../]
[]

[UserObjects]
  [./inserter]
    type = DiscreteNucleationInserter
    hold_time = 1e-2 #lower limit = 7.2e-5
    probability = P
    radius = r_crit
    time_dependent_statistics = False
  [../]

  [./mappsi]
    type = DiscreteNucleationMap
    periodic = cv
    inserter = inserter
    int_width = 1
  [../]

  [./term]
    type = Terminator
    expression = 'rate > 1'
  [../]
[]


[Postprocessors]
  [./time]
    type = TimePostprocessor
  [../]
  [./total_F]
    type = ElementIntegralVariablePostprocessor
    variable = f_dens
  [../]
  [./total_cv]
    type = ElementIntegralVariablePostprocessor
    variable = cv
    execute_on = 'initial timestep_end'
  [../]
  [./avg_c] #Need to calculate nucleation probability
    type = ElementAverageValue
    variable = cv
    execute_on = 'initial timestep_end'
  [../]
  [./diffusivity_avg] #Needed to calculate left BC
    type = ElementAverageMaterialProperty
    mat_prop = D
    execute_on = 'initial timestep_end'
  [../]
  [./length_scale] #Needed to calculate left BC
    type = ElementAverageMaterialProperty
    mat_prop = length_scale
    execute_on = 'initial timestep_end'
  [../]
  [./time_scale] #Needed to calculate left BC
    type = ElementAverageMaterialProperty
    mat_prop = time_scale
    execute_on = 'initial timestep_end'
  [../]
  [./Va] #Needed to calculate left BC
    type = ElementAverageMaterialProperty
    mat_prop = Va
    execute_on = 'initial timestep_end'
  [../]
  [./left_BC]
    type = FunctionValuePostprocessor
    function = Left_BC
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
    p2nucleus = 0.01
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
  nl_abs_tol = 1e-10
  l_tol = 1e-05
  nl_max_its = 15
  l_max_its = 30
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    cutback_factor = 0.8
    growth_factor = 1.2
    optimal_iterations = 7
  [../]

  [./Adaptivity]
    max_h_level = 2
    refine_fraction = 0.9
    coarsen_fraction = 0.05
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
