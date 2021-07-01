[StochasticTools]
[]

[Samplers]
  [hypercube]
    type = LatinHypercube
    num_rows = 8
    distributions = 'Em Ef flux_exponent gamma'
    num_bins = 5
  []
  [cartesian]
    type = CartesianProduct
    linear_space_items = '0.05 0.05 4
                          0.2  0.1  7
                          -5   1    6
                          0.1  0.6  4'
  []
  [csv]
    type = CSVSampler
    # samples_file = 'grid_out_minimal.csv'
    # samples_file = 'grid_out_fortest.csv'
    # samples_file = 'grid_out_stdcase.csv'
    # samples_file = 'grid_out_only_Em_vary_narrow_data_0000.csv'
    # samples_file = 'grid_out_only_Ef_vary_narrow_data_0000.csv'
    # samples_file = 'grid_out_only_flux_vary_narrow_data_0000.csv'
    samples_file = 'grid_out_only_sig_vary_narrow_data_0000.csv'
    execute_on = 'initial'
  []
[]

[GlobalParams]
  # sampler = hypercube
  # sampler = cartesian
  sampler = csv
[]

[Distributions]
  [Em]
    type = Uniform
    lower_bound = 0.055
    upper_bound = 0.2
  []
  [Ef]
    type = Uniform
    lower_bound = 0.055
    upper_bound = 0.55
  []
  [flux_exponent]
    type = Uniform
    lower_bound = -5
    upper_bound = 0
  []
  [gamma]
    type = Uniform
    lower_bound = 0.1
    upper_bound = 1.8
  []
[]

[MultiApps]
  [runner]
    type = SamplerFullSolveMultiApp
    #Defined in GlobalParams
    # sampler = hypercube
    # sampler = cartesian
    # sampler = csv
    input_files = 'heterogen_nucl_diff_only_with_sink_1D_ideal_sol.i'
    # mode = batch-reset
    mode = normal
    ignore_solve_not_converge = true
  []
[]

[Transfers]
  [parameters]
    type = SamplerParameterTransfer
    multi_app = runner
    #Defined in GlobalParams
    # sampler = hypercube
    # sampler = cartesian
    # sampler = csv
    parameters = 'Materials/Em/prop_values Materials/Ef/prop_values Materials/flux_exponent/prop_values Materials/gamma_Jm2/prop_values'
    to_control = 'stochastic'
  []
  [results]
    type = SamplerPostprocessorTransfer
    multi_app = runner
    #Defined in GlobalParams
    # sampler = hypercube
    # sampler = cartesian
    # sampler = csv
    to_vector_postprocessor = results
    from_postprocessor = 'time max_cv'
  []
[]

[VectorPostprocessors]
  [results]
    type = StochasticResults
  []
  [samples]
    type = SamplerData
    #Defined in GlobalParams
    # sampler = hypercube
    # sampler = cartesian
    # sampler = csv
  []
  # [stats]
  #   type = Statistics
  #   vectorpostprocessors = results
  #   compute = 'mean'
  #   ci_method = 'percentile'
  #   ci_levels = '0.05'
  # []
[]

[Outputs]
  csv = true
  # execute_on = 'FINAL'
  # file_base = 'master_out'
  # file_base = 'runtest_out'

  # file_base = 'Em_vary_out'
  # file_base = 'Ef_vary_out'
  # file_base = 'flux_vary_out'
  file_base = 'sig_vary_out'

[]
