[Mesh]
  type = GeneratedMesh
  xmax = 100 #nm
  # ymax = 100 #nm
  nx = 50
  # ny = 50
  dim = 1
[]

[Variables]
  [./cv]
  []
[]

[ICs]
  [./cv]
    type = ConstantIC
    variable = cv
    value = 2e-7
  [../]
[]

[Kernels]
  [./Diff]
    type = MatDiffusion
    variable = cv
    diffusivity = D
  [../]
  [./Cdot]
    type = TimeDerivative
    variable = cv
  [../]
[]

[Materials]
  [./Conditions] #Sets the temperature, length scale, and time scale of the simulation
    type = GenericConstantMaterial
    prop_names =  'T   length_scale time_scale'
    prop_values = '300 1e-9         1e-6      '
  [../]
  [./Material_properties] #Sets the free energy properties of the material
    type = GenericConstantMaterial
    prop_names =  'Em    D0'
    prop_values = '0.055 8.37707e-7 '
  [../]
  [./Constants]
    type = GenericConstantMaterial
    prop_names =  'kb'
    prop_values = '8.6173303e-5'
  [../]
  [./diffusivity]
    type = ParsedMaterial
    material_property_names = 'length_scale time_scale Em kb T D0'
    function = 'D0/length_scale^2*time_scale*exp(-Em/(kb*T))'
    f_name = D
    outputs = exodus
  [../]
  [./Va] #Atomic volume
    type = ParsedMaterial
    f_name = Va
    material_property_names = 'length_scale'
    constant_names = 'Va'
    constant_expressions = 20.24 #In angstroms^3/atom (Ma and Dudarev 2019)
    function = 'Va*(1e-10/length_scale)^3'  #Change units of length
  [../]
[]

[Functions]
  [./Left_BC]
    type = ParsedFunction
    #vars = flux (mol/m^2-s) length_scale time_scale Va (nm^3/atom) Na (atoms/mol)
    #The length_scale, time_scale and Va here should match in materials block
    vars = 'flux length_scale time_scale Va Na         D'
    vals = '1e-2 length_scale time_scale Va 6.02214e23 diffusivity_avg'
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
[]

[AuxVariables]
  [./analytical]
  [../]
[]

[AuxKernels]
  [./function]
    type = FunctionAux
    variable = analytical
    function = analytical
  [../]
[]

[BCs]
  [./left_flux]
    type = FunctionNeumannBC
    variable = cv
    boundary = left
    function = Left_BC
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_abs_tol = 1e-13
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    optimal_iterations = 3
    cutback_factor = 0.5
    growth_factor = 1.1
  [../]
[]

[VectorPostprocessors]
  [./cv_vec]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '100 0 0'
    sort_by = 'x'
    num_points = 100
    variable = cv
  [../]
  [./analytical_vec]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '100 0 0'
    sort_by = 'x'
    num_points = 100
    variable = analytical
  [../]
[]

[Postprocessors]
  [./left_value]
    type = SideAverageValue
    boundary = left
    variable = cv
  [../]
  [./right_value]
    type = SideAverageValue
    boundary = right
    variable = cv
  [../]
  [./dt]
    type = TimestepSize
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
  [./left_flux]
    type = SideFluxAverage
    variable = cv
    diffusivity = 1.0
    boundary = left
  [../]
  [./left_flux_an]
    type = SideFluxAverage
    variable = analytical
    diffusivity = 1.0
    boundary = left
  [../]
  [./left_val_an]
    type = SideAverageValue
    variable = analytical
    boundary = left
  [../]
  [./error]
    type = ElementL2Error
    function = analytical
    variable = cv
  [../]
  [./left_val_no_grad]
    type = FunctionValuePostprocessor
    function = no_gradient
  [../]
[]

[UserObjects]
  [./end]
    type = Terminator
    expression = 'left_value > 1e-3'
  [../]
[]

[Outputs]
  exodus = true
[]
