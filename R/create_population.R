




create_population_specific_ages <- function(N, initial_worms, contact_rates_by_age,
                                            worm_stages, female_factor, male_factor,initial_miracidia,
                                            initial_miracidia_days, predis_aggregation, time_step,
                                            spec_ages, ages_per_index, death_rate_per_time_step,
                                            mda_adherence, mda_access){
  
  args <- list(N, initial_worms, contact_rates_by_age,
               worm_stages, female_factor, male_factor,initial_miracidia,
               initial_miracidia_days, predis_aggregation, time_step,
               spec_ages, ages_per_index, death_rate_per_time_step,
               mda_adherence, mda_access)
  names(args) <- c('N', 'initial_worms', 'contact_rates_by_age',
                   'worm_stages', 'female_factor', 'male_factor', 'initial_miracidia',
                   'initial_miracidia_days', 'predis_aggregation', 'time_step',
                   'spec_ages', 'ages_per_index', 'death_rate_per_time_step',
                   'mda_adherence', 'mda_access') 
  
  JuliaCall::julia_assign("args", args)
  JuliaCall::julia_eval("result = create_population_specified_ages")
  
  
}