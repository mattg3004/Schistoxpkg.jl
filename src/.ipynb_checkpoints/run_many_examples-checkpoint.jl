using Schistoxpkg
using Random
using JLD
using Distributions
#using Plots

# include parameters and add a file name
 include("/Users/matthewgraham/Dropbox/Schistoxpkg/Schistoxpkg.jl/src/parameters.jl")
 filename = "equ.jld"

 # make parameters structure to hold all the parameters for the simulation
 pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
         density_dependent_fecundity, average_worm_lifespan,
         max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
         worm_stages, contact_rate, max_fec_contact_rate_product, max_fecundity, age_contact_rates,
         ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
         birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
         death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
         spec_ages, ages_per_index, record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
         rate_acquired_immunity, M0, human_larvae_maturity_time, egg_sample_size)
 pars = make_age_contact_rate_array(pars, scenario, [],[]);

min_max_fecundity = 50
max_max_fecundity = 500

min_max_fec_contact_rate = 0.5
max_max_fec_contact_rate = 2

min_predis_aggregation = 0.001
max_predis_aggregation = 0.08

min_sac_burden = 3
max_sac_burden = 15

#min_sac_heavy_burden = 38
#max_sac_heavy_burden = 48



min_cercariae_survival = 0.00001
max_cercariae_survival = 0.1

min_miracidia_survival = 0.00001
max_miracidia_survival = 0.1


mutable struct outRuns
    times
    prev
    sac_prev
    adult_prev
    high_burden
    high_burden_sac
    high_adult_burden
    max_fecundity
    max_fec_contact_rate_product
    predis_aggregation
    miracidia_survival
    cercariae_survival
end




function runAllSims(pars, min_max_fecundity, max_max_fecundity,
                min_max_fec_contact_rate, max_max_fec_contact_rate,
                min_predis_aggregation, max_predis_aggregation,
                min_miracidia_survival, max_miracidia_survival,
                min_sac_burden, max_sac_burden, number_years, numSims)

    allRunsRecord = Array{Union{Missing, outRuns}}(missing, numSims)

    for i in 1:numSims
        print("at run ", i)
        print("")
        output = run_1_sim(pars, min_max_fecundity, max_max_fecundity,
                    min_max_fec_contact_rate, max_max_fec_contact_rate,
                    min_predis_aggregation, max_predis_aggregation,
                    min_miracidia_survival, max_miracidia_survival,
                    min_sac_burden, max_sac_burden, number_years)

        allRunsRecord[i] =  output
    end
    return allRunsRecord
end



function runAllSimsParallel(pars, min_max_fecundity, max_max_fecundity,
                min_max_fec_contact_rate, max_max_fec_contact_rate,
                min_predis_aggregation, max_predis_aggregation,
                min_miracidia_survival, max_miracidia_survival,
                min_sac_burden, max_sac_burden, number_years, numSims)

    allRunsRecord = Array{Union{Missing, outRuns}}(missing, numSims)

    Threads.@threads for i in 1:numSims
        output = run_1_sim(pars, min_max_fecundity, max_max_fecundity,
                    min_max_fec_contact_rate, max_max_fec_contact_rate,
                    min_predis_aggregation, max_predis_aggregation,
                    min_miracidia_survival, max_miracidia_survival,
                    min_sac_burden, max_sac_burden, number_years)

        allRunsRecord[i] =  output
    end
    return allRunsRecord
end



function run_1_sim(pars, min_max_fecundity, max_max_fecundity,
                    min_max_fec_contact_rate, max_max_fec_contact_rate,
                    min_predis_aggregation, max_predis_aggregation,
                    min_miracidia_survival, max_miracidia_survival,
                    min_sac_burden, max_sac_burden, number_years)



        output = outRuns(0, 0, 0, 0, 0, 0, 0,
                                  0, 0,  0,
                                  0, 0)
                     # create the larvae variables along with the human structure
      humans, miracidia, cercariae = create_population_specified_ages(pars)

                     # update the ages and death ages, so that there aren't any people with the death age lower than actual age
      humans = generate_ages_and_deaths(20000, humans, pars)
      humans = update_contact_rate(humans,  pars)


                     #
      number_years = 100
      num_time_steps = trunc(Int, 365*number_years / time_step)

      pars.max_fecundity = min_max_fecundity + rand()*(max_max_fecundity-min_max_fecundity)
      pars.max_fec_contact_rate_product = min_max_fec_contact_rate + rand()*(max_max_fec_contact_rate-min_max_fec_contact_rate)
      pars.predis_aggregation = min_predis_aggregation + rand() *( max_predis_aggregation - min_predis_aggregation)

      pars.miracidia_survival = min_miracidia_survival + rand()*(max_miracidia_survival - min_miracidia_survival)
      pars.cercariae_survival = pars.miracidia_survival


      # cercariae_survival = runif(1, min_cercariae_survival, max_cercariae_survival)
      #update_parameters_individually(pars, name = "cercariae_survival", value = miracidia_survival)
      number_years_equ = 65
      times, sac_burden, sac_heavy_burden, record, humans, miracidia,cercariae = run_to_equ2(pars, number_years_equ)
      equ_sac_burden  = sac_burden[length(sac_burden)]
      equ_sac_heavy_burden  = sac_heavy_burden[length(sac_heavy_burden)]
      print("equ_sac_burden =", equ_sac_burden)
      print("equ_sac_heavy_burden =", equ_sac_heavy_burden)
      if (length(sac_burden) > 10)

          if ((equ_sac_burden > min_sac_burden) & (equ_sac_burden < max_sac_burden))

          num_repeats = 3 #number of simulations to run
          number_years = 20
          drug_efficacy = 0.863 #Toor et al. JID paper in SI: drug efficacy 86.3% for S. mansoni and 94% for S. haematobium
          num_time_steps =  trunc(Int, 365*number_years / time_step)

          mda_start1 = 1
          mda_end1 = 4.1

          mda_info = create_mda(0.75, .75, 0.75, mda_start1,
                                mda_end1, 1, [0,1], [0,1], [0,1], drug_efficacy)




          vaccine_info =  []
          times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden =
                run_with_MDA2(pars, humans, miracidia, cercariae, mda_info, num_time_steps, num_repeats)



                output = outRuns(times, prev, sac_prev, adult_prev, high_burden, high_burden_sac, high_adult_burden,
                              pars.max_fecundity,  pars.max_fec_contact_rate_product,  pars.predis_aggregation,
                              pars.miracidia_survival, pars.cercariae_survival)
        end

    end

    return output
end


function run_with_MDA2(pars, humans, miracidia, cercariae, mda_info,num_time_steps ,num_repeats )
  #  source("Initial_conditions.R")
  # contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)


  vaccine_info = []

  humans, miracidia, cercariae, record = update_env_constant_population(num_time_steps, humans,
                                                                              miracidia, cercariae, pars, mda_info, vaccine_info)




  sac_burden, high_burden_sac, times = get_sac_data_from_record(record)
  adult_prev, high_adult_burden, times = get_adult_data_from_record(record)
  prev, high_burden, times = get_population_data_from_record(record)


  return  times, prev, sac_burden, high_burden, high_burden_sac, adult_prev, high_adult_burden
end


function get_sac_data_from_record(record)
    sac_burden = (p->p.sac_burden[1]).(record)
    sac_heavy_burden = (p->p.sac_burden[3]).(record)
    times = (p->p.time).(record);
    return sac_burden, sac_heavy_burden, times
end


function get_adult_data_from_record(record)
    adult_burden = (p->p.adult_burden[1]).(record)
    adult_heavy_burden = (p->p.adult_burden[3]).(record)
    times = (p->p.time).(record);
    return adult_burden, adult_heavy_burden, times
end


function get_population_data_from_record(record)
      population_burden = (p->p.population_burden[1]).(record)
      population_heavy_burden = (p->p.population_burden[3]).(record)
      times = (p->p.time).(record);
      return population_burden, population_heavy_burden, times
end

function run_to_equ2(pars, number_years_equ)


  humans, miracidia, cercariae = create_population_specified_ages(pars)
  humans = generate_ages_and_deaths(20000, humans, pars)
  humans = update_contact_rate(humans,  pars)

  num_time_steps = trunc(Int, 365*number_years_equ / time_step)

  mda_info = []

  vaccine_info = []

   humans, miracidia, cercariae, record =
  @time update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info);

  # ages, eggs, female_worms, male_worms = get_ages_eggs_worms(humans)


  #multiSimsBinCounts = make_multiSimBins(ages, eggs, age_groups, bins, egg_multiplier)

  sac_burden, sac_heavy_burden, times = get_sac_data_from_record(record)

  return times, sac_burden, sac_heavy_burden, record, humans, miracidia,cercariae
end





Threads.threadid()
Threads.nthreads()
number_years = 65
numSims = 20

allRunsRecord = @time runAllSims(pars, min_max_fecundity, max_max_fecundity,
                min_max_fec_contact_rate, max_max_fec_contact_rate,
                min_predis_aggregation, max_predis_aggregation,
                min_miracidia_survival, max_miracidia_survival,
                min_sac_burden, max_sac_burden, number_years, numSims)


allRunsRecord = @time runAllSimsParallel(pars, min_max_fecundity, max_max_fecundity,
                                min_max_fec_contact_rate, max_max_fec_contact_rate,
                                min_predis_aggregation, max_predis_aggregation,
                                min_miracidia_survival, max_miracidia_survival,
                                min_sac_burden, max_sac_burden, number_years, numSims)
