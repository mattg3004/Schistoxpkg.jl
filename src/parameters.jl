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
worm_stages = 1
female_factor = 1
male_factor = 1
contact_rate = 0.03
ages_per_index = 5

# if more than one community, then specify how many here
#N_communities = 3
N_communities = 1

# next parameter is the relative probabilities of being in each community
# if entries are all equal, then all communities are equally likely and will
# be roughly the same size
#community_probs = [1,1,1]
community_probs = 1

# community_contact_rate = [1,1,1]
community_contact_rate = 1

# parameter for proportion of people who are given mda who will take it
mda_adherence = .9
mda_access = .9

heavy_burden_threshold = 16
# number of days after which miracidia become cercariae
miracidia_maturity_time = 24 # for S. mansoni
# miracidia_maturity_time = 21 # for S. haemotobium

env_cercariae = 0
initial_miracidia = 50000*N/1000
init_env_cercariae = 50000*N/1000
initial_miracidia_days = trunc(Int,ceil(miracidia_maturity_time/time_step, digits = 0))

# how long to run simulation for
number_years_equ = 200

max_fecundity = 0.34  # for S. mansoni [Toor et al JID paper SI]
#max_fecundity = 0.3  # for S. haematobium [Toor et al JID paper SI]

density_dependent_fecundity = 0.0007 # for S. mansoni [Toor et al JID paper SI]
#density_dependent_fecundity = 0.0006 # for S. haematobium [Toor et al JID paper SI]

r = 0.03 # aggregation parameter for negative binomial for egg production
num_time_steps_equ = trunc(Int, 365*number_years_equ / time_step)

# human birth rate
birth_rate = 28*time_step/(1000*365)

average_worm_lifespan = 5.7 # years for S. mansoni [Toor et al JID paper SI]
#average_worm_lifespan = 4 # years for S. haematobium [Toor et al JID paper SI]

# this is the aggregation parameter for the predisposition
predis_aggregation = 0.24
predis_weight = 1

# what proportion of miracidias and cercariae survive each round
env_miracidia_survival_prop = 1/2
env_cercariae_survival_prop = 1/2
mda_coverage = 0.8 # proportion of target age group reached by mda
mda_round = 0

# proportion of cercariae which can infect humans
human_cercariae_prop = 1

# gamma distribution for Kato-Katz method
gamma_k = Gamma(0.87,1/0.87)

vaccine_effectiveness = 0.95
num_sims = 1

# record the state of the population this often in years
record_frequency = 1/24

#= this is the number of thousands of people in 5 year (0-4, 5-9,...) intervals in Kenya
and will be used to give a specified age structure when we run to equilibrium =#
spec_ages = 7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
            2636, 1970, 1468, 1166, 943, 718, 455, 244

#= number of deaths per 1000 individuals by age
    first entry is for under 1's, then for 5 year intervals from then on =#


death_prob_by_age = [0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                                                0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1]

ages_for_deaths = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                                              65, 70, 75, 80, 85, 90, 95, 100, 110]


drug_efficacy = 0.863
scenario = "moderate adult"
