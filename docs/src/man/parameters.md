
# Parameters

```
N: size of human population (Int)

time_step: length of time in days to step forward each time (Float)

N_communities: number of communities sharing the same environmental source (Int)

community_probs: probability of an integer being in a given community. This governs
the population sizes in each community. This must be an array the length of
N_communities (Array(Float))

community_contact_rate: contact rate with the environment for each of the communities. 
This must be an array the length of N_communities (Array(Float))

density_dependent_fecundity: decrease in egg production per worm due to high density
of worms (Float)

average_worm_lifespan: average life expectancy of a worm (Float)

max_age: maximum age of human in the population (Float)

initial_worms: number of worms in each person to begin with. The actual number in 
each person if chosen from a Poisson distirubiton with this mean (Int)

initial_miracidia: initial number of miracidia larvae in the environment (Int)

initial_miracidia_days: miracidia will age into cercariae larvae after a specified
number of days. This parameter will specify how many days of initial_miracidia we 
will have already had in the environment (Int)

init_env_cercariae: initial number of cercariae larvae in the environment (Int)

worm_stages: how many age stages are there for the worms. Having 1 stage will give
Gamma distributed death ages, while more than 1 will result in Erlang distribution (Int)

contact_rate: global contact rate for the uptake of larvae from the environment (Float)

max_fec_contact_rate_product: product of maximum fecundity and contact rate parameters. This
is often used to keep the simulation in a desired range of behaviour when maximum fecundity is changed.
We will set the contact rate depending on this parameter and the maximum fecundity. Appropriate values
of this parameter will depend on the value of the cercariae and miracidia survival. (Float)

max_fecundity: expected number of eggs from a single worm pair. The actual number 
will be chosen from a distribution with this mean (Float)

age_contact_rates: contact rate for chosen age groups(Array(Float))

ages_for_contacts: age groups for specifying contact rates (Array(Int))

contact_rate_by_age_array: array holding contact rate for each age (Array(Float))

mda_adherence: proportion of people who adhere to the mda (Float)

mda_access: proportion of people who have access to the mda (Float)

female_factor: factor for altering the contact rate for females, if we choose to 
have gender specific behaviour which affects contact rate (Float)

male_factor: factor for altering the contact rate for males, if we choose to have 
gender specific behaviour which affects contact rate (Float)

miracidia_maturity: number of days after which miracidias will mature to cercariae (Int)

birth_rate: rate of birth of humans (Float)

human_cercariae_prop: proportion of cercariae which are able to infect humans (Float)

predis_aggregation: aggregation for predisposition of individuals to uptake larvae.
This is chosen from a Gamma distribution with mean 1 for each individual and set for 
life. If this is high, then the aggregation is low, meaning that most individuals
have roughly the same predisposition. If it is low, then the larvae become
concentrated in a few individuals. (Float)

cercariae_survival: what proportion of cercariae survive from one time point
to the next (Float)

miracidia_survival: what proportion of miracidia survive from one time point 
to the next (Float)

death_prob_by_age: for specified age range, what is the probability of dying 
each year (Array(Float))

ages_for_death: age ranges for death probabilities (Array(Float))

r: aggregation parameter for negative binomially distributed egg production (Float)

vaccine_effectiveness: efficacy of a vaccine if one is used(Float)

drug_effectiveness: efficacy of a drug given during MDA(Float)

spec_ages: number of individuals by age group which we specify if we want a 
particular age distribution for the simulation (Array(Float))

ages_per_index: how many different ages we include in the spec_ages parameter (Int)

record_frequency: how often we should record the prevalence in the population 
during simulation (Float)

use_kato_katz: if 0, then don't use Kato-Katz (KK) for egg counts, if 1, use KK (Int)

kato_katz_par: parameter for Gamma distribution if KK is used (Float)

heavy_burden_threshold: number of eggs at which an individual is said to have a
heavy infection (Int)

rate_acquired_immunity: rate at which immunity will be acquired for individuals.
This will be multiplied by the cumulative number of worms people have had throughout
their life to decide the level of immunity acquired (Float)

M0: if a particular for of egg production is used, this parameter is required and is
a proxy for mean worm burden (Float) 

human_larvae_maturity_time: length of time in days after which a cercariae uptaken by
a human will mature into a worm (Int)

egg_sample_size: the sample size of daily output for egg counting. For haematobium, this 
will be ~1/100 as the usual sample size is 10mL and a daily amount of urine is ~1L (Float >0 & <1)

egg_production_distribution: Either "NegBin" of "Poisson", depending on whether a negative 
binomial or Poisson distribution is desired for the production of eggs given the 
number of worms in an individual (String)

```
