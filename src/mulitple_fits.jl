contact_rate = 0.0051
r = 0.03 # aggregation parameter for negative binomial for egg production
human_cercariae_prop = 1
env_cercariae_survival_prop = 1/2
env_miracidia_survival_prop = 1/4
N = 2000
max_fecundity = 0.34
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [], []);
# contact_rates_by_age = make_age_contact_rate_array(max_age, scenario,[],[]);

# contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [4,9,15,20, 25, 100], [0.01, 0.5, 1.5, 0.5, 0.5, 0.01])
## run to equilibrium
ages , gender, predisposition,  human_cercariae, eggs, vac_status,
        treated, female_worms, male_worms, death_rate, age_contact_rate,
        vaccinated, env_miracidia, adherence, access = create_population_specified_ages(N, initial_worms, contact_rates_by_age,
                    worm_stages, female_factor, male_factor,initial_miracidia,
                    initial_miracidia_days, predis_aggregation, time_step,
                    spec_ages, ages_per_index, death_rate_per_time_step,
                    mda_adherence, mda_access)



ages_equ, gender_equ, predisposition_equ,  human_cercariae_equ, eggs_equ,
vac_status_equ, treated_equ, female_worms_equ, male_worms_equ,
vaccinated_equ, env_miracidia_equ, env_cercariae_equ, record_high =
    update_env_to_equilibrium(num_time_steps_equ, copy(ages), copy(human_cercariae), copy(female_worms), copy(male_worms),
        time_step, average_worm_lifespan,
        copy(eggs), max_fecundity, r, worm_stages,
        copy(vac_status), copy(gender), predis_aggregation,
        copy(predisposition), copy(treated), vaccine_effectiveness,
        density_dependent_fecundity, copy(vaccinated), copy(env_miracidia),
        copy(env_cercariae), contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
        female_factor, male_factor, contact_rates_by_age, record_frequency, copy(age_contact_rate),human_cercariae_prop)
