# Schistoxpkg.jl

A package to run an individual based mode of a schistosomiasis outbreak. 

All parameters are stored in the parameters.jl file in the src folder.

To run to equilibrium, a few functions need to be fun. After defining parameters we have to run "contact_rates_by_age" = make_age_contact_rate_array(max_age, scenario, [], []); to generate the contact rate by age array. After this, we have to run:

ages , death_ages, gender, predisposition, community,
   human_cercariae, eggs, vac_status,
   treated, female_worms, male_worms, age_contact_rate,
   vaccinated, env_miracidia, adherence, access = 
   create_population_specified_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                worm_stages, female_factor, male_factor,initial_miracidia,
                initial_miracidia_days, predis_aggregation, predis_weight,
                time_step,
                spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                mda_adherence, mda_access)

death_ages = []
for i in 1:N
     push!(death_ages, get_death_age(death_prob_by_age, ages_for_deaths))
end
mean(death_ages)
ages, death_ages = generate_ages_and_deaths(20000, ages, death_ages, death_prob_by_age, ages_for_deaths)

ages_equ, gender_equ, predisposition_equ,  human_cercariae_equ, eggs_equ,
vac_status_equ, treated_equ, female_worms_equ, male_worms_equ,
vaccinated_equ, env_miracidia_equ, env_cercariae_equ, record_high =
    update_env_to_equilibrium(7600, copy(ages), copy(human_cercariae), copy(female_worms), copy(male_worms),
    copy(community), community_contact_rate,
        time_step, average_worm_lifespan,
        copy(eggs), max_fecundity, r, worm_stages,
        copy(vac_status), copy(gender), predis_aggregation,
        copy(predisposition), copy(treated), vaccine_effectiveness,
        density_dependent_fecundity, copy(vaccinated), copy(env_miracidia),
        copy(env_cercariae), contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
        female_factor, male_factor, contact_rates_by_age, record_frequency, copy(age_contact_rate),human_cercariae_prop,
        miracidia_maturity_time)

```@index
```

```@autodocs
Modules = [Schistoxpkg]
```
