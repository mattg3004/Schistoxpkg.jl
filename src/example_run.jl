using Schistoxpkg
using Random
using JLD
using Distributions
using PyPlot

# include parameters and add a file name
include("parameters.jl")
filename = "equ.jld"

# make parameters structure to hold all the parameters for the simulation
pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
        worm_stages, contact_rate, max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
        death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index, record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold)
pars = make_age_contact_rate_array(pars, scenario, [],[]);

# create the larvae variables along with the human structure
humans, miracidia, cercariae = create_population_specified_ages(pars)

# update the ages and death ages, so that there aren't any people with the death age lower than actual age
humans = generate_ages_and_deaths(20000, humans, pars)
humans = update_contact_rate(humans,  pars)


#
number_years_mda = 200
num_time_steps = trunc(Int, 365*number_years_mda / time_step)
num_repeats = 1

# create mda information
# mda_info = create_mda(0, .75, 0, 1, number_years_mda, 1, [0,1], [0,1], [0,1], pars.drug_effectiveness)
mda_info = []

vaccine_info = []

humans, miracidia, cercariae, record =
        update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)
#@time update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info);

# save population
save_population_to_file(filename, humans,  miracidia, cercariae, pars)

# plot school age children prevalence and heavy burden
pp = (p->p.sac_burden[1]).(record)
pp2 = (p->p.sac_burden[3]).(record)
times = (p->p.time).(record);
clf()
plt.plot(times,pp)
plt.plot(times,pp2)
gcf()



# do simulations with mda for school age children every year for 20 years

num_repeats = 10 #number of simulations to run
number_years = 20
drug_efficacy = 0.863 #Toor et al. JID paper in SI: drug efficacy 86.3% for S. mansoni and 94% for S. haematobium
num_time_steps = 365*number_years / time_step

mda_info = create_mda(0, .75, 0, 1, number_years, 1, [0,1], [0,1], [0,1], pars.drug_effectiveness)


vaccine_info = []
times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden =
run_repeated_sims_no_births_deaths(filename, num_time_steps, mda_info, vaccine_info, num_repeats)
#run_repeated_sims_no_population_change(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

clf()
plt.plot(times, mean.(sac_prev))
plt.plot(times, mean.(high_burden_sac))
gcf()
