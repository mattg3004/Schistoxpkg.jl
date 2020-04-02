using JLD
using Plots
using Distributions
using Random

# filename to save population in =#
# filename = "equ_runs_1.jld"
N = 1000
max_age = 100
initial_worms = 10
time_step = 10
worm_stages = 2
female_factor = 1
male_factor = 1
initial_miracidia = 1
initial_miracidia_days = trunc(Int,round(41/time_step, digits = 0))
env_cercariae = 0
contact_rate = 0.00003
ages_per_index = 5

# parameter for proportion of people who are given mda who will take it
mda_adherence = .9
mda_access = .9

number_years = 250
max_fecundity = 0.34  # From "The design of schistosomiasis monitoring and evaluation programmes:
#The importance of collecting adu lt data to inform treatment strategies for Schistosoma mansoni"
density_dependent_fecundity = 0.0005
r = 0.03 # aggregation parameter for negative binomial for egg production
num_time_steps = trunc(Int, 365*number_years / time_step)

birth_rate = 28*time_step/(1000*365)

average_worm_lifespan = 5.7 # years

# this is the aggregation parameter for the predisposition
predis_aggregation = 0.24

env_cercariae_death_rate = 0.09 * time_step #= life span of cercariae in the environment is short 8-20 hrs
according to "Studies of the Transmission Dynamics, Mathematical Model Development and the Control of Schistosome
Parasites by Mass Drug Administration in Human Communities"  =#
env_miracidia_death_rate = 0 * time_step
mda_coverage = 0.8 # proportion of target age group reached by mda
mda_round = 0

# gamma distribution for Kato-Katz method
gamma_k = Gamma(0.87,1/0.87)

vaccine_effectiveness = 0.95
num_sims = 1

# record the state of the population this often in years
record_frequency = 1/12

#= this is the number of thousands of people in 5 year (0-4, 5-9,...) intervals in Kenya
and will be used to give a specified age structure when we run to equilibrium =#
spec_ages = 7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
            2636, 1970, 1468, 1166, 943, 718, 455, 244

#= number of deaths per 1000 individuals by age
    first entry is for under 1's, then for 5 year intervals from then on =#
age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                           21.83, 29.98, 36.98]
