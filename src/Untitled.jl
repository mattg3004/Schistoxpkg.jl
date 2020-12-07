

function multiple_runs_to_equ(total_number_runs, num_time_steps_equ, num_time_steps_mda, num_repeats_mda,
  pars, mda_info, vaccine_info, filename)

   times_equ = []
   all_prev_mda = []
   all_sac_prev_mda = []
   all_high_burden_mda = []
   all_high_burden_sac_mda =[]
   all_adult_prev_mda = []
   all_high_adult_burden_mda = []

   times_mda = []
   prev_equ = []
   sac_prev_equ = []
   high_burden_equ = []
   high_burden_sac_equ = []
   adult_prev_equ = []
   high_adult_burden_equ = []

   for run in 1 : total_number_runs

     humans, miracidia, cercariae = create_population_specified_ages(pars)

     # update the ages and death ages, so that there aren't any people with the death age lower than actual age
     humans = generate_ages_and_deaths(20000, humans, pars)
     humans = update_contact_rate(humans,  pars)


      humans, miracidia, cercariae, record = update_env_to_equilibrium(num_time_steps_equ, humans, miracidia, cercariae, pars)

      push!(prev_equ, (p->p.pop_prev).(record))
      push!(sac_prev_equ, (p->p.sac_burden[1]).(record))
      push!(high_burden_equ, (p->p.population_burden[3]).(record))
      push!(high_burden_sac_equ, (p->p.sac_burden[3]).(record))
      push!(adult_prev_equ, (p->p.adult_prev).(record))
      push!(high_adult_burden_equ, (p->p.adult_burden[3]).(record))
      times_equ  = (p ->p.time).(record)


      save_population_to_file(filename, humans,  miracidia, cercariae, pars)

      times_mda, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden =
              run_repeated_sims_no_births_deaths_human_larvae(filename, num_time_steps_mda, mda_info, vaccine_info, num_repeats_mda)

      push!(all_prev_mda, prev)
      push!(all_sac_prev_mda , sac_prev)
      push!(all_high_burden_mda , high_burden)
      push!(all_high_burden_sac_mda , high_burden_sac)
      push!(all_adult_prev_mda , adult_prev)
      push!(all_high_adult_burden_mda , high_adult_burden)

    end
      return prev_equ, sac_prev_equ, high_burden_equ, high_burden_sac_equ, adult_prev_equ, high_adult_burden_equ, times_equ,
      all_prev_mda, all_sac_prev_mda, all_high_burden_mda, all_high_burden_sac_mda, all_adult_prev_mda, all_high_adult_burden_mda, times_mda

   end

total_number_runs = 2
num_time_steps_equ = 3650
num_time_steps_mda = 730
num_repeats_mda = 20


prev_equ, sac_prev_equ, high_burden_equ, high_burden_sac_equ, adult_prev_equ, high_adult_burden_equ, times_equ,
all_prev_mda, all_sac_prev_mda, all_high_burden_mda, all_high_burden_sac_mda, all_adult_prev_mda, all_high_adult_burden_mda, times_mda =
    multiple_runs_to_equ(total_number_runs, num_time_steps_equ, num_time_steps_mda, num_repeats_mda,
      pars, mda_info, vaccine_info, filename)
