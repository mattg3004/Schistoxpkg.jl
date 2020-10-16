var documenterSearchIndex = {"docs":
[{"location":"#Schistoxpkg.jl","page":"Home","title":"Schistoxpkg.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package to run an individual based mode of a schistosomiasis outbreak based on original code from this paper. Generally people uptake larvae based on a contact rate defined by their age, along with some predisposition which is chosen from a gamma distribution with mean 1, but some specified level of variance.","category":"page"},{"location":"","page":"Home","title":"Home","text":"All parameters are stored in the parameters.jl file in the src folder.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The model has a parameter which defines the time step that we take forward each time. Several functions are then called each time step which simulate the run of the outbreak. This is repeated until we reach a specified number of steps, usually corresponding to stepping forward a chosen number of years into the future.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The standard approach is given by the following set of processes which all have their own function to execute them. First, load required packages and include the parameters file which stores the parameters for the model.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Schistoxpkg\nusing Random\nusing JLD\nusing Distributions\nusing PyPlot\ninclude(\"parameters.jl\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"After this, we must initialize the pars struct, which will be used to store the parameters.","category":"page"},{"location":"","page":"Home","title":"Home","text":"pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,\n                  density_dependent_fecundity, average_worm_lifespan,\n                  max_age, initial_worms, initial_miracidia, \n                  initial_miracidia_days, init_env_cercariae,\n                  worm_stages, contact_rate, max_fecundity, age_contact_rates,\n                  ages_for_contacts, contact_rate_by_age_array, mda_adherence, \n                  mda_access,  female_factor, male_factor, miracidia_maturity,\n                  birth_rate, human_cercariae_prop, predis_aggregation, \n                  cercariae_survival, miracidia_survival,\n                  death_prob_by_age, ages_for_death, r, \n                  vaccine_effectiveness, drug_effectiveness,\n                  spec_ages, ages_per_index, record_frequency, \n                  use_kato_katz, kato_katz_par, heavy_burden_threshold,\n                  rate_acquired_immunity, M0, human_larvae_maturity_time)\n\npars = make_age_contact_rate_array(pars, scenario, [],[]);","category":"page"},{"location":"","page":"Home","title":"Home","text":"We then counstruct the initial host population and initialize the environment cercariae and mircacidia larvae. We also age the population forward for a certain number of steps (here 20,000), so that there are no individuals in the population whose death age is lower than their actual age, and the death rate dynamics have time to achieve a desired age distribution. We also update the age based contact rates according to the new ages in the population.","category":"page"},{"location":"","page":"Home","title":"Home","text":"humans, miracidia, cercariae = create_population_specified_ages(pars)\nhumans = generate_ages_and_deaths(20000, humans, pars)\nhumans = update_contact_rate(humans,  pars)\n\nmda_info = []\n\nvaccine_info = []","category":"page"},{"location":"","page":"Home","title":"Home","text":"Each time step we advance the time of the simulation by the length of the time step and also add this time step to the age of each individual. There is a chosen period at which contact rates are updated for each individual, where we check if someone has aged into a different age bracket, resulting if their level of contact has changed.","category":"page"},{"location":"","page":"Home","title":"Home","text":"We then calculate the total number of worms within each individual and the number of pairs of worms a person has. These numbers are used to calculate how many eggs someone will produce. The number of eggs is chosen from a poisson distribution with mean equal to the number of worm pairs multiplied by the max fecundity parameter and then multiplied by an exponential function which calculates the density dependent reduction in eggs produced, λ wp exp(-wp z). We then kill the worms within human hosts at a given rate, which is based on average worm lifespan.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Eggs are then hatched into the environment, with egg release dependent on the age specific contact rate of each individual. Humans are given an age of death when they are born, which is based on some chosen death rates for each age group. We check each time step if anyone has outlived their age of death and if they have, they are then removed from the population. Cercariae from the environment are then uptaken to each surviving individual based on their predisposition and contact rate. These immediately become worms within the human host.","category":"page"},{"location":"","page":"Home","title":"Home","text":"We then perform any interventions which are due to take place at this point in time after which we will cull the miracidia and cercariae in the environment by a chosen percentage. After this we will add births to the population which occur at some specified rate.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A form of this is done by with the following function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"number_years = 200\nnum_time_steps = trunc(Int, 365*number_years / time_step)\nhumans, miracidia, cercariae, record = \nupdate_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, \n                                         cercariae, pars, mda_info, vaccine_info)","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are other versions of this basic approach, where we don't age the population or include births and deaths and also where the population is aged but every death is simply matched with a birth, resulting in the population being kept constant. dfhfdhaergewadsazrgtsw","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Schistoxpkg]","category":"page"},{"location":"#Schistoxpkg.cercariae_uptake!-NTuple{4,Any}","page":"Home","title":"Schistoxpkg.cercariae_uptake!","text":"cercariae_uptake(humans, cercariae, miracidia, pars)\n\nuptake cercariae into humans, whilst updating cercariae with miracidia. Uptaken cercariae immediately become worms in this formulation\n\n\n\n\n\n","category":"method"},{"location":"#Schistoxpkg.cercariae_uptake_with_human_larvae!-NTuple{4,Any}","page":"Home","title":"Schistoxpkg.cercariae_uptake_with_human_larvae!","text":"cercariae_uptake(humans, cercariae, miracidia, pars)\n\nuptake cercariae into humans, whilst updating cercariae with miracidia. Uptaken cercariae immediately become worms in this formulation\n\n\n\n\n\n","category":"method"},{"location":"#Schistoxpkg.create_population-Tuple{Any}","page":"Home","title":"Schistoxpkg.create_population","text":"create_population\n\nThis will create the initial human population with randomly chosen age, and gender. Predisposition is taken to be gamma distributed There is also a male and female adjustment to predisposition adjusting for gender specific behaviour In addition to this, it will create the initial miracidia environment vector\n\n\n\n\n\n","category":"method"},{"location":"#Schistoxpkg.create_population_specified_ages-Tuple{Any}","page":"Home","title":"Schistoxpkg.create_population_specified_ages","text":"create_population_specified_ages\n\nThis will create the initial human population with an age distribution specified by the spec_ages variable Predisposition is taken to be gamma distributed. There is also a male and female adjustment to predisposition adjusting for gender specific behaviour In addition to this, it will create the initial miracidia environment vector\n\n\n\n\n\n","category":"method"},{"location":"#Schistoxpkg.egg_production!-Tuple{Any,Any}","page":"Home","title":"Schistoxpkg.egg_production!","text":"egg_production(humans, pars)\n\nfunction to produce eggs for individuals, dependent on how many worms they have         and the max fecundity and density dependent fecundity of the population\n\n\n\n\n\n","category":"method"},{"location":"#Schistoxpkg.update_contact_rate-Tuple{Any,Any}","page":"Home","title":"Schistoxpkg.update_contact_rate","text":"update_contact_rate(humans, pars)\n\nfunction to update the contact rate of individuals in the population. This is necessary     as over time when people age, they will move through different age groups which have     different contact rates\n\n\n\n\n\n","category":"method"},{"location":"#Schistoxpkg.worm_maturity!-Tuple{Any,Any}","page":"Home","title":"Schistoxpkg.worm_maturity!","text":"worm_maturity(humans, pars)\n\nfunction to kill worms, and if there is more than one stage for worm life,         to update how many worms are in each stage\n\n\n\n\n\n","category":"method"}]
}
