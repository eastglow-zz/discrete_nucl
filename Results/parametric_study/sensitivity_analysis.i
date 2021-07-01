[StochasticTools]
[]

[Samplers]
  [hypercube]
    type = LatinHypercube
    num_rows = 8
    distributions = 'Em flux_exponent gamma'
    num_bins = 3
  []

  [hypercube_a]
    type = LatinHypercube
    # num_rows = 10000
    num_rows = 100
    distributions = 'Em flux_exponent gamma'
    num_bins = 3
    seed = 2011
  []
  [hypercube_b]
    type = LatinHypercube
    # num_rows = 10000
    num_rows = 100
    distributions = 'Em flux_exponent gamma'
    num_bins = 3
    seed = 2013
  []
  [sobol]
    type = Sobol
    sampler_a = hypercube_a
    sampler_b = hypercube_b
  []

[]

[GlobalParams]
  # sampler = hypercube
  # sampler = cartesian
  sampler = sobol
[]

[Distributions]
  [Em]
    type = Uniform
    lower_bound = 0.055
    upper_bound = 0.2
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
    input_files = 'heterogen_nucl_diff_only_with_sink.i'
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
    parameters = 'Materials/Em/prop_values Materials/flux_exponent/prop_values Materials/gamma_Jm2/prop_values'
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
  [sobol]
    type = SobolStatistics
    sampler = sobol
    results = results
  []
[]

[Outputs]
  csv = true
  # execute_on = 'FINAL'
  # file_base = 'master_out'
  # file_base = 'runtest_out'
  # file_base = 'sig_vary_out'
  # file_base = 'flux_vary_out'
  file_base = 'SA_out'
[]
