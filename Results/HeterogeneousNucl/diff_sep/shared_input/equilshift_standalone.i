#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 200 #Number of elements in the x-direction
  ny = 200 #Number of elements in the y-direction
  xmax = 100 #X-direction domain size, for this file it is in nm
  ymax = 100 #Y-direction domain size, for this file it is in nm
  # uniform_refine = 2
[]

[GlobalParams] #Parameters used in multiple input blocks
  fa_name  = fb
  fb_name  = fv
  ca       = cvB
  cb       = cvV
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

[AuxVariables]
  [./f_dens]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./psi_map]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # Identifies surface regions for heterogeneous nucleation
  [./surf]
  [../]
  # Identifies channel regions for heterogeneous nucleation
  [./channel]
  [../]
  # Storage variable for the average vacancy concentration
  [./cv_avg]
  [../]
  [./cvB_max]
  [../]
[]

[ICs]
  [./cv_IC]
    type = ConstantIC
    variable = cv
    value = 2.021e-7
    # value = 0.01
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
    # value = 0.01
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
    int_width = 0.6928
  [../]
  [./channel_IC]
    type = BoundingBoxIC
    variable = channel
    x1 = -1
    x2 = 101
    y1 = 48
    y2 = 52
    inside = 1
    outside = 0
    int_width = 0.6928
  [../]
[]

[Functions]
  [./get_max_cvB]
    type = ParsedFunction
    vars = 'max_cvB'
    vals = 'max_cvB'
    value = max_cvB
  [../]
  [./Left_BC]
    type = ParsedFunction
    #vars = flux (mol/m^2-s) length_scale time_scale Va (nm^2/atom) Na (atoms/mol)
    #The length_scale, time_scale and Va here should match in materials block
    vars = 'flux length_scale time_scale Va Na         D'
    vals = '1e-2 length_scale time_scale Va 6.02214e23 diffusivity_avg'
    # value = 'flux*Na*length_scale^2*time_scale*Va/D'
    value = 'flux*Na*length_scale^2*time_scale*Va'
  [../]
  [./Left_slope]
    type = ParsedFunction
    #vars = flux (mol/m^2-s) length_scale time_scale Va (nm^3/atom) Na (atoms/mol)
    #The length_scale, time_scale and Va here should match in materials block
    vars = 'flux length_scale time_scale Va Na         D'
    vals = '1e-2 length_scale time_scale Va 6.02214e23 diffusivity_avg'
    value = 'flux*Na*length_scale^2*time_scale*Va/D'
  [../]
  [./no_gradient]
    type = ParsedFunction
    vars = 'D dCdx C0 lx'
    vals = 'diffusivity_avg Left_slope 2e-7 100'
    value = 'C0 + dCdx*D/lx*t'
  [../]
  [./total_time_function]
    type = ParsedFunction
    vars = 'time_diff_only'
    vals = 'time_diff_only'
    value = 't + time_diff_only'
  [../]
[]

[Materials]
  ######Parameters defining the material behavior
  [./Conditions] #Sets the temperature, length scale, and time scale of the simulation
    type = GenericConstantMaterial
    prop_names =  'T   length_scale time_scale'  #Units: K, m, s
    prop_values = '300 1e-9         1e-6      '
  [../]
  [./Material_properties] #Sets the material properties of the material
    type = GenericConstantMaterial
    prop_names =  'Em    Ef   Ab    Bb eps  Av   cVeq D0          dcB   sink'
    prop_values = '0.055 0.52 270.5 32 0.01 1000 1    8.37707e-7  0.02     0'
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
    material_property_names = 'length_scale time_scale Em kb T'
    constant_names = 'D0'
    constant_expressions = '8.37707e-7'  # m^2/s
    # function = '1.39*(1e-2/length_scale^2)*time_scale*exp(-Em/(kb*T))'
    function = 'D0/length_scale^2*time_scale*exp(-Em/(kb*T))'
    f_name = D
    outputs = exodus
  [../]
  [./Mobility]
    type = ParsedMaterial
    f_name = M
    material_property_names = 'D h fVcc:=D[fv,cvV,cvV] fBcc:=D[fb,cvB,cvB]'
    function = 'fcc:=fVcc*fBcc/((1-h)*fVcc+h*fBcc + 1e-30);D/fcc'
    outputs = exodus
  [../]
  [./PF_mob]
    type = ParsedMaterial
    f_name = L
    material_property_names = 'L0 D kappa'
    function = 'L0*D/kappa'
    outputs = exodus
  [../]

  ####### Calculation of free energy expressions
  [./fb_fit]
    type = DerivativeParsedMaterial
    args = cvB
    material_property_names = 'Ab Bb eps cBeq'
    function = 'Ab*(sqrt((cvB-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(cvB-cBeq)^2'
    f_name = fb
    outputs = exodus
  [../]
  [./fv_poly]
    type = DerivativeParsedMaterial
    args = cvV
    material_property_names = 'Av cVeq'
    function = '1/2 * Av * (cvV - cVeq)^2'
    f_name = fv
    outputs = exodus
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
    material_property_names = 'delta gamma'
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
    # penalty = 160
    penalty_mode = MIN
    map = mappsi
    outputs = exodus
  [../]

  [./prob_homeg]
    type = ParsedMaterial
    f_name = P_homog
    args = 'cv psi'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em P_filter'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(P_filter<0.01,0,P_filter)*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    outputs = exodus
  [../]

  [./prob_heterog]
    type = ParsedMaterial
    f_name = P_heterog
    args = 'cv psi'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em S P_filter'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(P_filter<0.01,0,P_filter)*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
    # function = '0.0001'
    outputs = exodus
  [../]

  [./probability]
    type = ParsedMaterial
    f_name = P
    args = 'cv psi surf'
    material_property_names = 'P_homog P_heterog'
    # function = '(surf*P_heterog + (1-surf)*P_homog)'
    function = 'P_homog'
    outputs = exodus
  [../]

  [./shape_factor] #Heterogeneous nucleation shape factor
    type = ParsedMaterial
    f_name = S
    material_property_names = 'pi'
    constant_names =       'theta'
    constant_expressions = '60'
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
    function = 'if(dgv<0,-1,1)*pi*gamma^2/(abs(dgv)+1e-30)'  #in 2D
    # function = '16*pi*gamma^3/3/(abs(dgv))^2'  #in 3D
    f_name = dgstar
    outputs = exodus
  [../]

  [./rstar] #Critical radius
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
    function = 'Ab*(sqrt((cv-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(cv-cBeq)^2'
    args = 'cv'
    f_name = dgvmatrix
    outputs = exodus
  [../]

  [./dgvnucl]
    type = ParsedMaterial
    material_property_names = 'Av cVeq'
    args = 'cv_avg'
    function = '1/2 * Av * (cVeq - cVeq)^2' #This sets the energy at nucleation, and currently is set to zero
    f_name = dgvnucl
    outputs = exodus
  [../]

  [./fcbulk_slope]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb eps'
    args = 'cv'
    function = 'Ab*(cv/sqrt((cv-cBeq)^2+eps^2)) + Bb*cv'
    f_name = fcbulk_slope
    outputs = exodus
  [../]

  [./fcbulk_tangent]
    type = ParsedMaterial
    material_property_names = 'cVeq dgvmatrix fcbulk_slope'
    args = 'cv'
    function = 'fcbulk_slope*(cVeq-cv) + dgvmatrix'
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
    function = 'exp(Sf)*exp(-Ef/kb/T)' #From Frank et al. PRL 1996
    f_name = cBeq
    outputs = exodus
  [../]

  [./nucl_radius]
    type = ParsedMaterial
    function = '1.5'
    f_name = r_crit
  [../]

  [./cvB_max]
    type = ParsedMaterial
    args = 'psi cv'
    function = 'if(psi>1e-7,0,1)*cv'
    f_name = cvB_max_mat
    outputs = exodus
  [../]

  [./vac_sink_mask]
    type = ParsedMaterial
    args = 'channel cv'
    material_property_names = 'sink'
    function = 'if(cv<0,0,channel*sink)'
    f_name = vac_sink_mask
    outputs = exodus
  [../]

  [./segregation_strength]
    type = DerivativeParsedMaterial
    args = 'cv'
    material_property_names = 'cBeq dcB'
    constant_names = 'Bseg wd'
    constant_expressions = '1e3 0.05'
    function = '-Bseg*(0.5+0.5*tanh((cBeq+dcB - (cv+3*wd))/wd))'
    f_name = A_seg
    outputs = exodus
  [../]

  [./P_filter]
    type = ParsedMaterial
    args = 'cv psi'
    material_property_names = 'cBeq'
    f_name = P_filter
    function = 'if(psi>1e-7,0,1)*if(cv<0,0, cv)'
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
  [./sink]
    type = MaskedBodyForce
    variable = cv
    mask = vac_sink_mask
  [../]

  [./segregation]
    type = MatDiffusion
    variable = cv
    v = channel
    diffusivity = A_seg
    args = 'cv'
  [../]

  #
  # Allen-Cahn Equation
  #
  [./ACBulkF]
    type = KKSACBulkF
    variable = psi
    args     = 'cv cvB cvV'
    w        = 88.27 #kappa/delta^2
    mob_name = L
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = psi
    args = 'cv cvB cvV'
    mob_name = L
  [../]
  [./ACInterface]
    type = ACInterface
    variable = psi
    kappa_name = kappa
    mob_name = L
  [../]
  [./detadt]
    type = TimeDerivative
    variable = psi
  [../]

  [./fnuclpsi]
    type = AllenCahn
    variable = psi
    mob_name = L
    f_name = fnpsi
  [../]
  # [./forcenuclpsi]
  #   type = DiscreteNucleationForce
  #   variable = psi
  #   map = mappsi
  #   nucleus_value = 1
  #   no_nucleus_value = 0
  #   # mob_name = L
  #   # f_name = fnpsi
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
  [../]

  # [./av_cv]
  #   type = FunctionAux
  #   variable = cv_avg
  #   function = get_av
  # [../]

  [./max_cvB]
    type = FunctionAux
    variable = cvB_max
    function = get_max_cvB
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
    hold_time = 5 #lower limit = 7.2e-5
    probability = P
    radius = r_crit
    #seed = 230329348
    seed = 809982342 #Seed can be changed to change the random nucleation behavior
  [../]
  [./mappsi]
    type = DiscreteNucleationMap
    periodic = psi
    inserter = inserter
    int_width = 0.5
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
  [./avg_c] #Needed to calculate nucleation probability
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

  [./max_cvB]
    type = ElementExtremeMaterialProperty
    mat_prop = cvB_max_mat
    value_type = max
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
    p2nucleus = 0.01 #Increasing this will increase the time step before nucleation
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

  [./time_diff_only]
    type = Receiver
    # execute_on = INITIAL
  [../]

  [./total_time]
    type = FunctionValuePostprocessor
    function = total_time_function
    execute_on = 'initial timestep_end'
  [../]

  [./pfilter]
    type = ElementExtremeMaterialProperty
    value_type = min
    mat_prop = P_filter
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
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # petsc_options_iname = '-pc_type -ksp_type'
  # petsc_options_value = 'bjacobi  gmres'

  scheme = bdf2
  end_time = 10000

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  l_tol = 1e-04
  nl_max_its = 20
  l_max_its = 30
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    cutback_factor = 0.5
    growth_factor = 1.2
    optimal_iterations = 7
    linear_iteration_ratio = 100
    timestep_limiting_postprocessor = dtnuc
  [../]

  # [./Adaptivity]
  #   max_h_level = 2
  #   refine_fraction = 0.9
  #   coarsen_fraction = 0.05
  # [../]
[]

# [Adaptivity]
#   [./Indicators]
#     [./defect_channel]
#       type = GradientJumpIndicator
#       variable = channel
#     [../]
#     [./psigrad_jump]
#       type = GradientJumpIndicator
#       variable = psi
#     [../]
#     [./cvgrad_jump]
#       type = GradientJumpIndicator
#       variable = cv
#     [../]
#     [./cvBgrad_jump]
#       type = GradientJumpIndicator
#       variable = cvB
#     [../]
#     [./cvVgrad_jump]
#       type = GradientJumpIndicator
#       variable = cvV
#     [../]
#   [../]
#   [./Markers]
#     [./nuc]
#       type = DiscreteNucleationMarker
#       map = mappsi
#     [../]
#     [./gradpsi]
#       type = ValueThresholdMarker
#       variable = psigrad_jump
#       coarsen = 0.1
#       refine = 0.5
#     [../]
#     [./gradcv]
#       type = ValueThresholdMarker
#       variable = cvgrad_jump
#       coarsen = 0.1
#       refine = 0.5
#     [../]
#     [./gradcvB]
#       type = ValueThresholdMarker
#       variable = cvBgrad_jump
#       coarsen = 0.1
#       refine = 0.5
#     [../]
#     [./gradcvV]
#       type = ValueThresholdMarker
#       variable = cvVgrad_jump
#       coarsen = 0.1
#       refine = 0.5
#     [../]
#     [./gradchannel]
#       type = ValueThresholdMarker
#       variable = defect_channel
#       coarsen = 0.1
#       refine = 0.5
#     [../]
#     [./combo]
#       type = ComboMarker
#       markers = 'nuc gradcv gradcvB gradcvV gradchannel'
#     [../]
#   [../]
#   max_h_level = 2
#   marker = combo
# []

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
  file_base = 'equil_shift_combined2_out'
[]

[MultiApps]
  [./Initial_diffusion]
    type = FullSolveMultiApp
    app_type = combinedApp
    execute_on = 'INITIAL'
    positions = '0 0 0'
    input_files = ./equilshift_diff_only.i
  [../]
[]

[Transfers]
  [./Get_time_InitDiff_to_PF]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = Initial_diffusion
    from_postprocessor = time
    to_postprocessor = time_diff_only
    execute_on = SAME_AS_MULTIAPP
    reduction_type = minimum
  [../]
  [./Get_cv_InitDiff_to_PF]
    type = MultiAppMeshFunctionTransfer
    direction = from_multiapp
    multi_app = Initial_diffusion
    # source_variable = cv_avg
    source_variable = cv
    variable = cv
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Get_cvV_InitDiff_to_PF]
    type = MultiAppMeshFunctionTransfer
    direction = from_multiapp
    multi_app = Initial_diffusion
    source_variable = cvV
    variable = cvV
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Get_cvB_InitDiff_to_PF]
    type = MultiAppMeshFunctionTransfer
    direction = from_multiapp
    multi_app = Initial_diffusion
    # source_variable = cv_avg
    source_variable = cv
    variable = cvB
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Get_UPDATE_InitDiff_to_PF]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = Initial_diffusion
    from_postprocessor = update
    to_postprocessor = update
    execute_on = SAME_AS_MULTIAPP
    reduction_type = minimum
  [../]
  [./Get_COUNT_InitDiff_to_PF]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = Initial_diffusion
    from_postprocessor = count
    to_postprocessor = count
    execute_on = SAME_AS_MULTIAPP
    reduction_type = minimum
  [../]
  [./Pass_psi_PF_to_InitDiff]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    multi_app = Initial_diffusion
    source_variable = psi
    variable = psi
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_cv_PF_to_InitDiff]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    multi_app = Initial_diffusion
    source_variable = cv
    variable = cv
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_cvV_PF_to_InitDiff]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    multi_app = Initial_diffusion
    source_variable = cvV
    variable = cvV
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./Pass_cvB_PF_to_InitDiff]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    multi_app = Initial_diffusion
    source_variable = cvB
    variable = cvB
    execute_on = SAME_AS_MULTIAPP
  [../]
[]
