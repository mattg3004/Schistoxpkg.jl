# Schistoxpkg.jl

A package to run an individual based model of a schistosomiasis outbreak based on original code from this [paper](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-019-3749-4). Generally people uptake larvae based on a contact rate defined by their age, along with
some predisposition which is chosen from a gamma distribution with mean 1, but some specified level of variance.

All parameters are stored in the parameters.jl file in the src folder.

The model has a parameter which defines the time step that we take forward each time. Several functions are then called each
time step which simulate the run of the outbreak. This is repeated until we reach a specified number of steps, usually corresponding
to stepping forward a chosen number of years into the future.




The standard approach is given by the following set of processes which all have their own function to execute them. First, load required packages and include the parameters file which stores the parameters for the model.
```
using Schistoxpkg
using Random
using JLD
using Distributions
using PyPlot
include("parameters.jl")
```
After this, we must initialize the pars struct, which will be used to store the parameters.

```
pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
                  density_dependent_fecundity, average_worm_lifespan,
                  max_age, initial_worms, initial_miracidia, 
                  initial_miracidia_days, init_env_cercariae,
                  worm_stages, contact_rate, max_fecundity, age_contact_rates,
                  ages_for_contacts, contact_rate_by_age_array, mda_adherence, 
                  mda_access,  female_factor, male_factor, miracidia_maturity,
                  birth_rate, human_cercariae_prop, predis_aggregation, 
                  cercariae_survival, miracidia_survival,
                  death_prob_by_age, ages_for_death, r, 
                  vaccine_effectiveness, drug_effectiveness,
                  spec_ages, ages_per_index, record_frequency, 
                  use_kato_katz, kato_katz_par, heavy_burden_threshold,
                  rate_acquired_immunity, M0, human_larvae_maturity_time)

pars = make_age_contact_rate_array(pars, scenario, [],[]);
```

We then counstruct the initial host population and initialize the environment cercariae and mircacidia larvae. We also age the population forward for a certain number of steps (here 20,000), so that there are no individuals in the population whose death age is lower than their actual age, and the death rate dynamics have time to achieve a desired age distribution. We also update the age based contact rates according to the new ages in the population.
```
humans, miracidia, cercariae = create_population_specified_ages(pars)
humans = generate_ages_and_deaths(20000, humans, pars)
humans = update_contact_rate(humans,  pars)

mda_info = []

vaccine_info = []
```

Each time step we advance the time of the simulation by the length of the time step and also add this time step to the age of each individual.
There is a chosen period at which contact rates are updated for each individual, where we check if someone has aged into a different age bracket, resulting if their
level of contact has changed.

We then calculate the total number of worms within each individual and the number of pairs of worms a person has.
These numbers are used to calculate how many eggs someone will produce. The number of eggs is chosen from a poisson distribution with mean equal to the
number of worm pairs multiplied by the max fecundity parameter and then multiplied by an exponential function which calculates the density dependent reduction in eggs produced,
`λ wp exp(-wp z)`.
We then kill the worms within human hosts at a given rate, which is based on average worm lifespan.

Eggs are then hatched into the environment, with egg release dependent on the age specific contact rate of each individual.
Humans are given an age of death when they are born, which is based on some chosen death rates for each age group. We check each time step if anyone has outlived their age of death and if they have, they are then removed from the population.
Cercariae from the environment are then uptaken to each surviving individual based on their predisposition and contact rate. These immediately become worms within the human host.

We then perform any interventions which are due to take place at this point in time after which we will cull the miracidia and cercariae in the environment by a chosen percentage. After this we will add births to the population which occur at some specified rate.

A version of this is done by with the following function:
```
number_years = 200
num_time_steps = trunc(Int, 365*number_years / time_step)
humans, miracidia, cercariae, record = 
update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, 
                                         cercariae, pars, mda_info, vaccine_info)
```

There are other versions of this basic approach, where we don't age the population or include births and deaths and also where the population is aged but every death is simply matched with a birth, resulting in the population being kept constant.

## Functions

```@index
```

```@autodocs
Modules = [Schistoxpkg]
```
