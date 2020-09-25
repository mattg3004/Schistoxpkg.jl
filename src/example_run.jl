using Schistoxpkg
using Random
using JLD
using Distributions
using PyPlot

# include parameters and add a file name
include("parameters.jl")
filename = "equ.jld"

# make contact rates by age array
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [], []);

# create population
ages , death_ages, gender, predisposition, community,
   human_cercariae, eggs, vac_status,
   treated, female_worms, male_worms, age_contact_rate,
   vaccinated, env_miracidia, adherence, access = create_population_specified_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                worm_stages, female_factor, male_factor,initial_miracidia,
                initial_miracidia_days, predis_aggregation, predis_weight,
                time_step,
                spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                mda_adherence, mda_access)

# add death ages for each person
death_ages = []
for i in 1:N
    push!(death_ages, get_death_age(death_prob_by_age, ages_for_deaths))
end

# make it so that no ones age is greater than their deatha age
ages, death_ages = generate_ages_and_deaths(20000, ages, death_ages, death_prob_by_age, ages_for_deaths)

# update contact rates
age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)

# run to equilibrium for time length defined in parameters.jl
ages_equ, gender_equ, predisposition_equ,  human_cercariae_equ, eggs_equ,
vac_status_equ, treated_equ, female_worms_equ, male_worms_equ,
vaccinated_equ, env_miracidia_equ, env_cercariae_equ, record_high =
    update_env_to_equilibrium(num_time_steps_equ, copy(ages), copy(human_cercariae), copy(female_worms), copy(male_worms),
    copy(community), community_contact_rate,
        time_step, average_worm_lifespan,
        copy(eggs), max_fecundity, r, worm_stages,
        copy(vac_status), copy(gender), predis_aggregation,
        copy(predisposition), copy(treated), vaccine_effectiveness,
        density_dependent_fecundity, copy(vaccinated), copy(env_miracidia),
        copy(env_cercariae), contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
        female_factor, male_factor, contact_rates_by_age, record_frequency, copy(age_contact_rate),human_cercariae_prop,
        miracidia_maturity_time)


# save variables to file
save(filename, "ages", ages_equ ,  "gender", gender_equ,"predisposition",   predisposition_equ,
            "human_cercariae", human_cercariae_equ, "eggs", eggs_equ, "community", community,
            "vac_status", vac_status_equ,"treated", treated_equ, "female_worms",  female_worms_equ, "male_worms", male_worms_equ,
            "vaccinated", vaccinated_equ,  "age_contact_rate", age_contact_rate, "death_ages", death_ages,
            "env_miracidia",env_miracidia_equ, "env_cercariae", env_cercariae_equ, "adherence", adherence, "access", access)


###### run baseline case where mda isn't missed ########

# decide how long to run mda experiments for, along with how many times to run for
number_years_mda = 20
num_time_steps = trunc(Int, 365*number_years_mda / time_step)
num_repeats = 10

# create mda information
mda_info = create_mda(0, .75, 0, 1, number_years_mda, 1, [0,1], [0,1], [0,1], drug_efficacy)

vaccine_info = []


# run simulations for mda analysis
times, prev, sac_prev, high_burden, high_burden_sac, adult_prev = run_repeated_sims_no_population_change(num_repeats, num_time_steps,
            time_step, average_worm_lifespan,community_contact_rate, community_probs,
            max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
            density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
            female_factor, male_factor, contact_rates_by_age,
            death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
            record_frequency, filename, human_cercariae_prop, miracidia_maturity_time)

# generate variables to plot
mean_prev_baseline = mean.(prev)
mean_sac_prev_baseline = mean.(sac_prev)
mean_high_burden_baseline = mean.(high_burden)
mean_high_burden_sac_baseline = mean.(high_burden_sac)
mean_adult_prev_baseline = mean.(adult_prev)

# clear figure and plot
clf()
plt.plot(times, mean_sac_prev_baseline)
gcf()
