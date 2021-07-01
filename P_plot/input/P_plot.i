#
# KKS void problem in the split form
#
[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100 #Number of elements in the x-direction
  ny = 100 #Number of elements in the y-direction
  xmin = 0.05  #Minimum Em  (eV)
  xmax = 0.06 #Maximum Em  (eV)
  ymin = 0.4  #Minimum sig  (J/m^2)
  ymax = 0.5 #Maximum sig  (J/m^2)
[]

[GlobalParams] #Parameters used in multiple input blocks
  fa_name  = fb
  fb_name  = fv
  #ca       = cvB
  #cb       = cvV
[]

[AuxVariables]
  # vacancy concentration
  [./cv]
    # initial_condition = 2.021e-7
    initial_condition = 1e-3
  [../]

  # bulk phase concentration (matrix)
  [./cvB]
    initial_condition = 2.021e-7
  [../]

  # hydrogen phase concentration (delta phase)
  [./cvV]
    initial_condition = 1
  [../]

  # [./x]
  # [../]
  # [./y]
  # [../]

  [./P_aux]
  [../]
[]

[ICs]
[../]

[Functions]
  [./prob_fnc]
    type = ParsedFunction
    vars = 'Em    sig   Va P_expdgstar P_prefactor kb T eVpJ length_scale sig_J2ev'
    # vals = '0.055 0.406 Va P_expdgstar P_prefactor kb T eVpJ length_scale sig_J2ev'
    vals = 'xc    0.406  yc  P_expdgstar P_prefactor kb T eVpJ length_scale sig_J2ev'
    value = 'P_prefactor/Va*exp(-Em/kb/T)*pow(P_expdgstar, (sig*sig_J2ev)^2)'
  [../]

  [./xc]
    type = ParsedFunction
    value = 'x'
  [../]
  [./yc]
    type = ParsedFunction
    value = 'y'
  [../]
  [./tc]
    type = ParsedFunction
    value = 't'
  [../]
[]

[AuxKernels]
  [./prob_fnc_auxK]
    type = FunctionAux
    variable = P_aux
    function = prob_fnc
  [../]
  [./cv_as_t_axis]
    type = ParsedAux
    variable = cv
    use_xyzt = true
    function = 't'
  [../]
  [./cvB_copied_from_cv]
    type = ParsedAux
    variable = cvB
    args = cv
    function = 'cv'
  [../]
  # [./x_auxK]
  #   type = ParsedAux
  #   use_xyzt = true
  #   variable = x
  #   function = 'x'
  # [../]
  # [./y_auxK]
  #   type = ParsedAux
  #   use_xyzt = true
  #   variable = y
  #   function = 'y'
  # [../]
[]

[Materials]

  #######Parameters defining the material behavior
  [./Conditions] #Sets the temperature, length scale, and time scale of the simulation
    type = GenericConstantMaterial
    prop_names =  'T   length_scale time_scale'
    prop_values = '300 1e-9         1         '
  [../]
  [./Material_properties] #Sets the free energy properties of the material
    type = GenericConstantMaterial
    prop_names =  'Em    Ef   Ab    Bb eps  Av   cVeq D0        '
    prop_values = '0.055 0.52 270.5 32 0.01 1000 1    8.37707e-7'
    outputs = 'exodus'
  [../]
  [./Model_parameters]
    type = GenericConstantMaterial
    prop_names =  'L0   delta'
    prop_values = '0.1  0.6928'
  [../]
  [./Constants]
    type = GenericConstantMaterial
    prop_names =  'eVpJ           kb           h           pi'
    prop_values = '6.24150934e+18 8.6173303e-5 4.1357e-15  3.141592'
  [../]

  [./Va] #Atomic volume
    type = ParsedMaterial
    f_name = Va
    material_property_names = 'length_scale'
    constant_names = 'Va_in'
    constant_expressions = 20.24 #In angstroms^3/atom (Ma and Dudarev 2019)
    function = 'Va_in*(1e-10/length_scale)^3'  #Change units of length
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

  [./sig_J2ev]
    type = ParsedMaterial
    material_property_names = 'length_scale eVpJ T '
    function = 'eVpJ*length_scale^2'
    f_name = sig_J2ev
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
    material_property_names = 'D h fVcc:=D[fv,cvV,cvV] fBcc:=D[fb,cvB,cvB]'
    function = 'fcc:=fVcc*fBcc/((1-h)*fVcc+h*fBcc + 1e-30);D/fcc'
    outputs = exodus
  [../]

  ####### Calculation of free energy expressions
  [./fb_fit]
    type = DerivativeParsedMaterial
    args = 'cvB'
    material_property_names = 'Ab Bb eps cBeq'
    function = 'Ab*(sqrt((cvB-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(cvB-cBeq)^2'
    # material_property_names = 'Av cBeq'
    # function = '1/2 * Av * (cvB - cBeq)^2'
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

  [./prob_homeg]
    type = ParsedMaterial
    f_name = P_homog
    args = 'cv'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    # function = 'if(psi>0.01,0,1)*if(cv_analytical<0,0,1)*cv_analytical*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    function = 'cv*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    outputs = exodus
  [../]

  [./P_expdgstar]
    type = ParsedMaterial
    f_name = P_expdgstar
    material_property_names = 'kb T dgv pi'
    function = 'exp(-pi/(abs(dgv)+1e-30)/kb/T)'
    outputs = exodus
  [../]

  [./P_prefactor]
    type = ParsedMaterial
    f_name = P_prefactor
    args = 'cv'
    material_property_names = 'Ns length_scale time_scale T'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    # function = 'if(psi>0.01,0,1)*if(cv_analytical<0,0,1)*cv_analytical*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    function = 'cv*kb_h*T*time_scale'
    outputs = exodus
  [../]

  [./probability]
    type = ParsedMaterial
    f_name = P
    material_property_names = 'P_homog'
    function = 'P_homog'
    outputs = exodus
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
    function = 'Ab*(sqrt((cv-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(cv-cBeq)^2'
    args = 'cv'
    # material_property_names = 'Av cBeq'
    # args = 'cv'
    # function = '1/2 * Av * (if(cv<0,0,cv) - cBeq)^2'
    f_name = dgvmatrix
    outputs = exodus
  [../]

  [./dgvnucl]
    type = ParsedMaterial
    material_property_names = 'Av cVeq'
    args = 'cv'
    function = '1/2 * Av * (cVeq - cVeq)^2' #This sets the energy at nucleation, and currently is set to zero
    f_name = dgvnucl
    outputs = exodus
  [../]

  [./fcbulk_slope]
    type = ParsedMaterial
    material_property_names = 'cBeq Ab Bb eps'
    args = 'cv'
    function = 'Ab*(cv/sqrt((cv-cBeq)^2+eps^2)) + Bb*cv'
    # material_property_names = 'Av cBeq'
    # args = 'cv'
    # function = 'Av * (if(cv<0,0,cv) - cBeq)'
    f_name = fcbulk_slope
    outputs = exodus
  [../]

  [./fcbulk_tangent]
    type = ParsedMaterial
    material_property_names = 'cVeq dgvmatrix fcbulk_slope'
    args = 'cv'
    function = 'fcbulk_slope*(cVeq-cv) + dgvmatrix'
    # args = cv
    # function = 'fcbulk_slope*(cVeq-if(cv<0,0,cv)) + dgvmatrix'
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
[]

[UserObjects]
[]

[Executioner]
  # type = Steady
  type = Transient
  #cv value evolution along t-axis
  start_time = 2e-7
  end_time = 0.1
  dt = 1e-4
[]

[Postprocessors]
  [./length_scale]
    type = ElementAverageMaterialProperty
    mat_prop = length_scale
    execute_on = 'initial timestep_end'
  [../]
  [./time_scale]
    type = ElementAverageMaterialProperty
    mat_prop = time_scale
    execute_on = 'initial timestep_end'
  [../]
  [./Va]
    type = ElementAverageMaterialProperty
    mat_prop = Va
    execute_on = 'initial timestep_end'
  [../]
  [./kb]
    type = ElementAverageMaterialProperty
    mat_prop = kb
    execute_on = 'initial timestep_end'
  [../]
  [./T]
    type = ElementAverageMaterialProperty
    mat_prop = T
    execute_on = 'initial timestep_end'
  [../]
  [./h]
    type = ElementAverageMaterialProperty
    mat_prop = h
    execute_on = 'initial timestep_end'
  [../]
  [./pi]
    type = ElementAverageMaterialProperty
    mat_prop = pi
    execute_on = 'initial timestep_end'
  [../]
  [./eVpJ]
    type = ElementAverageMaterialProperty
    mat_prop = eVpJ
    execute_on = 'initial timestep_end'
  [../]
  [./dgv]
    type = ElementAverageMaterialProperty
    mat_prop = dgv
    execute_on = 'initial timestep_end'
  [../]
  [./Ceq]
    type = ElementAverageMaterialProperty
    mat_prop = cBeq
    execute_on = 'initial timestep_end'
  [../]
  [./P_expdgstar]
    type = ElementAverageMaterialProperty
    mat_prop = P_expdgstar
    execute_on = 'initial timestep_end'
  [../]
  [./P_prefactor]
    type = ElementAverageMaterialProperty
    mat_prop = P_prefactor
    execute_on = 'initial timestep_end'
  [../]
  [./sig_J2ev]
    type = ElementAverageMaterialProperty
    mat_prop = sig_J2ev
    execute_on = 'initial timestep_end'
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
[]
