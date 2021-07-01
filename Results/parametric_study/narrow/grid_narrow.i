[StochasticTools]
[]

[Samplers/sample]
  type = CartesianProduct
  # linear_space_items = '0.05 0.05 4
  #                       -5   1    6
  #                       0.1  0.57  4'
  # linear_space_items = '0.05 0.075 3
  #                       -5   2.5   3
  #                       0.1  0.85  3'  #minimal set
  linear_space_items = '0.055 0.01 1
                        0.5   0.1  1
                        -2   1     1
                        0.1 0.05  7'  #
  #Standard case: Em = 0.055
  #               Ef = 0.5
  #               flux = -2
  #               sig = 0.25
  #full space: 0.55 0.01 16
  #            0.2  0.1  7
  #            -5   1    6
  #            0.1  0.05 7
  execute_on = 'initial timestep_end'
[]

[VectorPostprocessors/data]
  type = SamplerData
  sampler = sample
  execute_on = 'initial timestep_end'
[]

[Outputs]
  execute_on = 'INITIAL'
  csv = true
  # file_base = 'grid_out_only_Em_vary_narrow'
  # file_base = 'grid_out_only_Ef_vary_narrow'
  # file_base = 'grid_out_only_flux_vary_narrow'
  file_base = 'grid_out_only_sig_vary_narrow'
[]
