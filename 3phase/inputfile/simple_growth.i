#
# KKS void problem in the split form
#
[Mesh]
  type = GeneratedMesh
  dim = 2

  #with adaptivity
  nx = 50 #Number of elements in the x-direction
  ny = 50 #Number of elements in the y-direction
  xmax = 100 #X-direction domain size, for this file it is in nm
  ymax = 100 #Y-direction domain size, for this file it is in nm
  uniform_refine = 2

  #no adaptivity
  # nx = 200 #Number of elements in the x-direction
  # ny = 200 #Number of elements in the y-direction
  # xmax = 100 #X-direction domain size, for this file it is in nm
  # ymax = 100 #Y-direction domain size, for this file it is in nm

[]

[GlobalParams] #Parameters used in multiple input blocks
  # fa_name  = fb
  # fb_name  = fv
  # ca       = cvB
  # cb       = cvV
  # wi        = 88.2714
[]

[Variables]
  # order parameter for bulk phase
  [./psiB]
  [../]
  # order parameter for void phase
  [./psiV]
  [../]
  # order parameter for defect channel
  [./psiD]
  [../]
  # vacancy concentration
  [./cv]
  [../]
  # bulk phase concentration (matrix)
  [./cvB]
  [../]
  # hydrogen phase concentration (void phase)
  [./cvV]
  [../]
  # hydrogen phase concentration (defect channel)
  [./cvD]
  [../]
  # Lagrange multiplier
  [./lambda]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
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
[]

[ICs]
  [./cv_IC]
    type = SmoothCircleIC
    variable = cv
    x1 = 50
    y1 = 50
    invalue = 1
    outvalue = 0.05
    int_width = 0.6928
    radius = 5
  [../]
  # [./cv_IC]
  #   variable = cv
  #   type = SmoothCircleIC
  #   x1 = 21.0
  #   y1 = 73.0
  #   radius = 2
  #   invalue = 1
  #   outvalue = 2.021e-7
  #   int_width = 0.6928
  # [../]
  [./cvB_IC]
    type = SmoothCircleIC
    variable = cvB
    x1 = 50
    y1 = 50
    invalue = 2.021e-7
    outvalue = 2.021e-7
    int_width = 0.6928
    radius = 5
  [../]
  [./cvD_IC]
    type = ConstantIC
    variable = cvD
    value = 0.05
    # value = 0.05
  [../]
  [./cvV_IC]
    type = SmoothCircleIC
    variable = cvV
    x1 = 50
    y1 = 50
    invalue = 1
    outvalue = 1
    int_width = 0.6928
    radius = 5
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
  [./psiV_IC]
    type = SmoothCircleIC
    variable = psiV
    x1 = 50
    y1 = 50
    invalue = 1
    outvalue = 0
    int_width = 0.6928
    radius = 5
  [../]
  [./psiD_IC]
    type = SmoothCircleIC
    variable = psiD
    x1 = 50
    y1 = 50
    invalue = 1
    outvalue = 1
    int_width = 0.6928
    radius = 5
  [../]
  [./psiB_IC]
    type = SmoothCircleIC
    variable = psiB
    x1 = 50
    y1 = 50
    invalue = 0
    outvalue = 1
    int_width = 0.6928
    radius = 5
  [../]

  # [./psiB_IC]
  #   variable = psiB
  #   type = SmoothCircleIC
  #   x1 = 21.0
  #   y1 = 73.0
  #   radius = 2
  #   invalue = 0
  #   outvalue = 1
  #   int_width = 0.6928
  # [../]
  # [./psiV_IC]
  #   variable = psiV
  #   type = SmoothCircleIC
  #   x1 = 21.0
  #   y1 = 73.0
  #   radius = 2
  #   invalue = 1
  #   outvalue = 0
  #   int_width = 0.6928
  # [../]
[]

[Functions]
  [./get_av]
    type = ParsedFunction
    vars = 'avg_cvB'
    vals = 'avg_cvB'
    value = avg_cvB
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
    prop_names =  'Em    Ef   Ab    Bb eps  Av    cVeq D0          dcB   sink'
    prop_values = '0.055 0.52 270.5 32 0.01 1000  1    8.37707e-7  0.05  0'
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
    # material_property_names = 'Ab Bb eps cBeq'
    # function = 'Ab*(sqrt((cvB-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(cvB-cBeq)^2'
    material_property_names = 'Av cBeq'
    function = '1/2 * Av * (cvB - cBeq)^2'
    f_name = fb
    outputs = exodus
  [../]
  [./fd]
    type = DerivativeParsedMaterial
    args = cvD
    # material_property_names = 'Ab Bb eps cBeq dcB'
    # function = 'Ab*(sqrt((cvD-(cBeq+dcB))^2 + eps^2) - eps) + 0.5*Bb*(cvD-(cBeq+dcB))^2'
    material_property_names = 'Av cBeq dcB'
    function = '1/2 * Av * (cvD - cBeq - dcB)^2'
    f_name = fd
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
  # [./h]
  #   type = SwitchingFunctionMaterial
  #   function_name = h
  #   h_order = HIGH
  #   eta = psi
  # [../]
  # Switching functions for each phase
  # h1(eta1, eta2, eta3)
  [./hB]
    type = SwitchingFunction3PhaseMaterial
    eta_i = psiB
    eta_j = psiV
    eta_k = psiD
    f_name = hB
  [../]
  # h2(eta1, eta2, eta3)
  [./hD]
    type = SwitchingFunction3PhaseMaterial
    eta_i = psiD
    eta_j = psiV
    eta_k = psiB
    f_name = hD
  [../]
  # h3(eta1, eta2, eta3)
  [./hV]
    type = SwitchingFunction3PhaseMaterial
    eta_i = psiV
    eta_j = psiB
    eta_k = psiD
    f_name = hV
  [../]

  # Coefficients for diffusion equation
  [./DhB]
    type = DerivativeParsedMaterial
    material_property_names = 'D hB'
    function = 'D*hB'
    f_name = DhB
  [../]
  [./DhD]
    type = DerivativeParsedMaterial
    material_property_names = 'D hD'
    function = 'D*hD'
    f_name = DhD
  [../]
  [./DhV]
    type = DerivativeParsedMaterial
    material_property_names = 'D hV'
    function = 'D*hV'
    f_name = DhV
  [../]

  # [./g]
  #   type = BarrierFunctionMaterial
  #   function_name = g
  #   g_order = SIMPLE
  #   eta = psi
  # [../]

  # Barrier functions for each phase
  [./gB]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = psiB
    function_name = gB
  [../]
  [./gD]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = psiD
    function_name = gD
  [../]
  [./gV]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = psiV
    function_name = gV
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
    op_names  = psiB
    op_values = 0.0
    penalty = 88.2714
    penalty_mode = MIN
    map = mappsi
    outputs = exodus
  [../]

  [./prob_homeg]
    type = ParsedMaterial
    f_name = P_homog
    args = 'cv psiV'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(psiV>0.01,0,1)*if(cv>0,cv,0)*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    # function = 'if(cv>0,cv,0)*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(dgstar)/kb/T)'
    outputs = exodus
  [../]

  [./prob_heterog]
    type = ParsedMaterial
    f_name = P_heterog
    args = 'cv psiV'
    material_property_names = 'Ns dgstar length_scale time_scale kb T  Em S'
    constant_names =       'kb_h'
    constant_expressions = '2.084e10'
    function = 'if(psiV>0.01,0,1)*if(cv>0,cv,0)*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
    # function = 'if(cv>0,cv,0)*Ns*kb_h*T*time_scale*exp(-Em/(kb*T))*exp(-abs(S*dgstar)/kb/T)'
    # function = '0.0001'
    outputs = exodus
  [../]

  [./probability]
    type = ParsedMaterial
    f_name = P
    args = 'cv psiV surf'
    material_property_names = 'P_homog P_heterog'
    # function = '(surf*P_heterog + (1-surf)*P_homog)'
    function = 'P_homog'  #for test
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
    function = 'pi*gamma^2/(abs(dgv)+1e-30)'  #in 2D
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
    # material_property_names = 'cBeq Ab Bb eps'
    # function = 'Ab*(sqrt((if(cv_avg<0,0,cv_avg)-cBeq)^2 + eps^2) - eps) + 0.5*Bb*(if(cv_avg<0,0,cv_avg)-cBeq)^2'
    # args = 'cv_avg'
    material_property_names = 'Av cBeq'
    args = 'cv'
    function = '1/2 * Av * (if(cv<0,0,cv) - cBeq)^2'
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
    # material_property_names = 'cBeq Ab Bb eps'
    # args = 'cv_avg'
    # function = 'Ab*(if(cv_avg<0,0,cv_avg)/sqrt((if(cv_avg<0,0,cv_avg)-cBeq)^2+eps^2)) + Bb*if(cv_avg<0,0,cv_avg)'
    material_property_names = 'Av cBeq'
    args = 'cv'
    function = 'Av * (if(cv<0,0,cv) - cBeq)'
    f_name = fcbulk_slope
    outputs = exodus
  [../]

  [./fcbulk_tangent]
    type = ParsedMaterial
    material_property_names = 'cVeq dgvmatrix fcbulk_slope'
    args = 'cv'
    function = 'fcbulk_slope*(cVeq-if(cv<0,0,cv)) + dgvmatrix'
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

  [./cvB_avg]
    type = ParsedMaterial
    args = 'psiV cv'
    function = 'if(psiV>0.01,0,1)*cv'
    f_name = cvB_avg
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

[]

[Kernels]

  # 3-phase CH
  [./diff_time]
    type = TimeDerivative
    variable = cv
  [../]
  [./diff_cvB]
    type = MatDiffusion
    variable = cv
    diffusivity = DhB
    v = cvB
    args = 'psiB psiD psiV'
  [../]
  [./diff_cvD]
    type = MatDiffusion
    variable = cv
    diffusivity = DhD
    v = cvD
    args = 'psiB psiD psiV'
  [../]
  [./diff_cvV]
    type = MatDiffusion
    variable = cv
    diffusivity = DhV
    v = cvV
    args = 'psiB psiD psiV'
  [../]



  #
  # Allen-Cahn Equation  #Ongoing
  #
  # [./ACBulkF]
  #   type = KKSACBulkF
  #   variable = psi
  #   args     = 'cv cvB cvV'
  #   w        = 88.27 #kappa/delta^2
  #   mob_name = L
  # [../]
  # [./ACBulkC]
  #   type = KKSACBulkC
  #   variable = psi
  #   args = 'cv cvB cvV'
  #   mob_name = L
  # [../]
  # [./ACInterface]
  #   type = ACInterface
  #   variable = psi
  #   kappa_name = kappa
  #   mob_name = L
  # [../]
  # [./detadt]
  #   type = TimeDerivative
  #   variable = psi
  # [../]
  # Kernels for Allen-Cahn equation for psiB
  [./dpsiBdt]
    type = TimeDerivative
    variable = psiB
  [../]
  [./ACBulkFB]
    type = KKSMultiACBulkF
    variable  = psiB
    Fj_names  = 'fb fd fv'
    hj_names  = 'hB hD hV'
    gi_name   = gB
    eta_i     = psiB
    wi        = 88.2714
    args      = 'cvB cvD cvV psiD psiV'
  [../]
  [./ACBulkCB]
    type = KKSMultiACBulkC
    variable  = psiB
    Fj_names  = 'fb fd fv'
    hj_names  = 'hB hD hV'
    cj_names  = 'cvB cvD cvV'
    eta_i     = psiB
    args      = 'psiD psiV'
  [../]
  [./ACInterfaceB]
    type = ACInterface
    variable = psiB
    kappa_name = kappa
  [../]
  [./multiplerB]
    type = MatReaction
    variable = psiB
    v = lambda
    mob_name = L
  [../]

  # Kernels for Allen-Cahn equation for psiD
  [./dpsiDdt]
    type = TimeDerivative
    variable = psiD
  [../]
  # [./ACBulkFD]
  #   type = KKSMultiACBulkF
  #   variable  = psiD
  #   Fj_names  = 'fb fd fv'
  #   hj_names  = 'hB hD hV'
  #   gi_name   = gD
  #   eta_i     = psiD
  #   wi        = 88.2714
  #   args      = 'cvB cvD cvV psiB psiV'
  # [../]
  # [./ACBulkCD]
  #   type = KKSMultiACBulkC
  #   variable  = psiD
  #   Fj_names  = 'fb fd fv'
  #   hj_names  = 'hB hD hV'
  #   cj_names  = 'cvB cvD cvV'
  #   eta_i     = psiD
  #   args      = 'psiB psiV'
  # [../]
  # [./ACInterfaceD]
  #   type = ACInterface
  #   variable = psiD
  #   kappa_name = kappa
  # [../]
  # [./multiplerD]
  #   type = MatReaction
  #   variable = psiD
  #   v = lambda
  #   mob_name = L
  # [../]

  # Kernels for the Lagrange multiplier equation
  [./mult_lambda]
    type = MatReaction
    variable = lambda
    mob_name = 3
  [../]
  [./mult_ACBulkF_B]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    gi_name   = gB
    eta_i     = psiB
    wi        = 88.2714
    mob_name  = 1
    args      = 'cvB cvD cvV psiD psiV'
  [../]
  [./mult_ACBulkC_B]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    cj_names  = 'cvB cvD cvV'
    eta_i     = psiB
    args      = 'psiD psiV'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_B]
    type = SimpleCoupledACInterface
    variable = lambda
    v = psiB
    kappa_name = kappa
    mob_name = 1
  [../]
  [./mult_ACBulkF_D]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    gi_name   = gD
    eta_i     = psiD
    wi        = 88.2714
    mob_name  = 1
    args      = 'cvB cvD cvV psiB psiV'
  [../]
  [./mult_ACBulkC_D]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    cj_names  = 'cvB cvD cvV'
    eta_i     = psiD
    args      = 'psiB psiV'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_D]
    type = SimpleCoupledACInterface
    variable = lambda
    v = psiD
    kappa_name = kappa
    mob_name = 1
  [../]
  [./mult_ACBulkF_V]
    type = KKSMultiACBulkF
    variable  = lambda
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    gi_name   = gV
    eta_i     = psiV
    wi        = 88.2714
    mob_name  = 1
    args      = 'cvB cvD cvV psiB psiD'
  [../]
  [./mult_ACBulkC_V]
    type = KKSMultiACBulkC
    variable  = lambda
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    cj_names  = 'cvB cvD cvV'
    eta_i     = psiV
    args      = 'psiB psiD'
    mob_name  = 1
  [../]
  [./mult_CoupledACint_V]
    type = SimpleCoupledACInterface
    variable = lambda
    v = psiV
    kappa_name = kappa
    mob_name = 1
  [../]

  # # Kernels for Allen-Cahn equation for psiD
  # [./dpsiDdt]
  #   type = TimeDerivative
  #   variable = psiD
  # [../]

  # Kernels for constraint equation psiB + psiV + psiD = 1
  # eta3 is the nonlinear variable for the constraint equation
  [./psiVreaction]
    type = MatReaction
    variable = psiV
    mob_name = 1
  [../]
  [./psiBreaction]
    type = MatReaction
    variable = psiV
    v = psiB
    mob_name = 1
  [../]
  # [./psiDreaction]
  #   type = MatReaction
  #   variable = psiV
  #   v = psiD
  #   mob_name = 1
  # [../]
  [./one]
    type = BodyForce
    variable = psiV
    value = -1.0
  [../]

  # Phase concentration constraints
  [./chempotBD]
    type = KKSPhaseChemicalPotential
    variable = cvB
    cb       = cvD
    fa_name  = fb
    fb_name  = fd
  [../]
  [./chempotDV]
    type = KKSPhaseChemicalPotential
    variable = cvD
    cb       = cvV
    fa_name  = fd
    fb_name  = fv
  [../]
  # [./chempotVB]
  #   type = KKSPhaseChemicalPotential
  #   variable = cvV
  #   cb       = cvB
  #   fa_name  = fv
  #   fb_name  = fb
  # [../]
  # [./phaseconcentration]
  #   type = KKSMultiPhaseConcentration
  #   variable = cvV
  #   cj = 'cvB cvD cvV'
  #   hj_names = 'hB hD hV'
  #   etas = 'psiB psiD psiV'
  #   c = cv
  # [../]
  [./phaseconcentration]
    type = KKSMultiPhaseConcentration
    variable = cvV
    cj = 'cvB cvD cvV'
    hj_names = 'hB hD hV'
    etas = 'psiB psiD psiV'
    c = cv
  [../]

  #
  # [./fnuclpsi]
  #   type = AllenCahn
  #   variable = psiB
  #   mob_name = L
  #   f_name = fnpsi
  # [../]

  # [./psiB_nucl_force]
  #   type = DiscreteNucleationForce
  #   variable = psiB
  #   map = mappsi
  #   no_nucleus_value = 1
  #   nucleus_value = 0
  # [../]
  # [./psiV_nucl_force]
  #   type = DiscreteNucleationForce
  #   variable = psiV
  #   map = mappsi
  #   no_nucleus_value = 0
  #   nucleus_value = 1
  # [../]

[]

[AuxKernels]
  # [./f_dens]
  #   variable = f_dens
  #   type = KKSGlobalFreeEnergy
  #   w = 88.27 #kappa/delta^2
  # [../]
  [./Energy_total]
    type = KKSMultiFreeEnergy
    Fj_names  = 'fb  fd  fv'
    hj_names  = 'hB  hD  hV'
    gj_names  = 'gB  gD  gV'
    variable = f_dens
    w = 88.27
    interfacial_vars =  'psiB  psiD  psiV'
    kappa_names =       'kappa kappa kappa'
  [../]

  [./psi_nucl_map]
    type = DiscreteNucleationAux
    map = mappsi
    variable = psi_map
  [../]

  [./av_cv]
    type = FunctionAux
    variable = cv_avg
    function = get_av
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
    periodic = psiB
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
  [./avg_cvB]
    type = ElementAverageMaterialProperty
    mat_prop = cvB_avg
    execute_on = 'initial timestep_end'
  [../]
  [./porosity]
    type = ElementAverageValue
    variable = psiV
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

  # petsc_options_iname = '-pc_type -ksp_type'   #doesn't work
  # petsc_options_value = 'bjacobi  gmres'

  # petsc_options_iname = '-pc_type -sub_pc_type   -sub_pc_factor_shift_type'
  # petsc_options_value = 'asm       ilu            nonzero'

  # petsc_options_iname = '-pc_type -sub_pc_type'
  # petsc_options_value = 'asm lu'

  # scheme = bdf2
  end_time = 10000
  # num_steps = 10

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

[Adaptivity]
  [./Indicators]
    [./cvgrad_jump]
      type = GradientJumpIndicator
      variable = cv
    [../]
    [./cvBgrad_jump]
      type = GradientJumpIndicator
      variable = cvB
    [../]
    [./cvDgrad_jump]
      type = GradientJumpIndicator
      variable = cvD
    [../]
    [./cvVgrad_jump]
      type = GradientJumpIndicator
      variable = cvV
    [../]
    [./psiDgrad_jump]
      type = GradientJumpIndicator
      variable = psiD
    [../]
  [../]
  [./Markers]
    [./nuc]
      type = DiscreteNucleationMarker
      map = mappsi
    [../]
    [./gradcv]
      type = ValueThresholdMarker
      variable = cvgrad_jump
      coarsen = 0.05
      refine = 0.9
    [../]
    [./gradcvB]
      type = ValueThresholdMarker
      variable = cvBgrad_jump
      coarsen = 0.05
      refine = 0.9
    [../]
    [./gradcvD]
      type = ValueThresholdMarker
      variable = cvDgrad_jump
      coarsen = 0.05
      refine = 0.9
    [../]
    [./gradcvV]
      type = ValueThresholdMarker
      variable = cvVgrad_jump
      coarsen = 0.05
      refine = 0.9
    [../]
    [./gradpsiD]
      type = ValueThresholdMarker
      variable = psiDgrad_jump
      coarsen = 0.2
      refine = 0.5
    [../]

    [./combo]
      type = ComboMarker
      markers = 'nuc gradcv gradcvB gradcvD gradcvV gradpsiD'
    [../]
  [../]
  max_h_level = 2
  marker = combo
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
  file_base = 'growth_in_defect_out'
[]

# [MultiApps]
#   [./Initial_diffusion]
#     type = FullSolveMultiApp
#     app_type = combinedApp
#     execute_on = 'INITIAL'
#     positions = '0 0 0'
#     input_files = ./heterogen_nucl_diff_only_optimized.i
#   [../]
# []
#
# [Transfers]
#   [./Get_time_InitDiff_to_PF]
#     type = MultiAppPostprocessorTransfer
#     direction = from_multiapp
#     multi_app = Initial_diffusion
#     from_postprocessor = time
#     to_postprocessor = time_diff_only
#     execute_on = SAME_AS_MULTIAPP
#     reduction_type = minimum
#   [../]
#   [./Get_cv_InitDiff_to_PF]
#     type = MultiAppMeshFunctionTransfer
#     direction = from_multiapp
#     multi_app = Initial_diffusion
#     # source_variable = cv_avg
#     source_variable = cv_analytical
#     variable = cv
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Get_cvV_InitDiff_to_PF]
#     type = MultiAppMeshFunctionTransfer
#     direction = from_multiapp
#     multi_app = Initial_diffusion
#     source_variable = cvV
#     variable = cvV
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Get_cvB_InitDiff_to_PF]
#     type = MultiAppMeshFunctionTransfer
#     direction = from_multiapp
#     multi_app = Initial_diffusion
#     # source_variable = cv_avg
#     source_variable = cv_analytical
#     variable = cvB
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Get_UPDATE_InitDiff_to_PF]
#     type = MultiAppPostprocessorTransfer
#     direction = from_multiapp
#     multi_app = Initial_diffusion
#     from_postprocessor = update
#     to_postprocessor = update
#     execute_on = SAME_AS_MULTIAPP
#     reduction_type = minimum
#   [../]
#   [./Get_COUNT_InitDiff_to_PF]
#     type = MultiAppPostprocessorTransfer
#     direction = from_multiapp
#     multi_app = Initial_diffusion
#     from_postprocessor = count
#     to_postprocessor = count
#     execute_on = SAME_AS_MULTIAPP
#     reduction_type = minimum
#   [../]
#   [./Pass_psi_PF_to_InitDiff]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = Initial_diffusion
#     source_variable = psiV
#     variable = psi
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_cv_PF_to_InitDiff]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = Initial_diffusion
#     source_variable = cv
#     variable = cv
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_cvV_PF_to_InitDiff]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = Initial_diffusion
#     source_variable = cvV
#     variable = cvV
#     execute_on = SAME_AS_MULTIAPP
#   [../]
#   [./Pass_cvB_PF_to_InitDiff]
#     type = MultiAppMeshFunctionTransfer
#     direction = to_multiapp
#     multi_app = Initial_diffusion
#     source_variable = cvB
#     variable = cvB
#     execute_on = SAME_AS_MULTIAPP
#   [../]
# []
