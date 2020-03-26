module Schistoxpkg
include("parameters.jl")
include("functions.jl")

export create_population
export make_age_contact_rate_array
export make_death_rate_array
export birth_of_human
export cercariae_uptake
export worm_maturity
export calculate_worm_pairs
export calculate_total_worms
export egg_production
export miracidia_production
export death_of_human
export birth_of_human
export human_cercariae_maturity
export update_env
export find_death_rate
export cercariae_death
export miracidia_death
export mda
export administer_drug
export update_contact_rate
export update_death_rate
export update_vaccine
export vaccinate
export update_mda
export load_population_from_file
export kato_katz
export get_prevalences
export save_population_to_file
export run_simulation_from_loaded_population
export find_mean_worms_by_age
export run_simulation
export generate_age_distribution
export specified_age_distribution
export create_population_specified_ages
export update_env_to_equilibrium
export create_mda
export vaccine_information
export mda_information
export update_env_no_births_deaths
export update_env_keep_population_same
export collect_prevs
export run_repeated_sims_no_population_change
export run_repeated_sims_random_births_deaths
export create_contact_settings


# using Distributions
# using Random
# using PoissonRandom

end # module
