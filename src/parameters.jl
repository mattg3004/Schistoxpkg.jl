using JLD
using Distributions
using Random

N = 1000    #population size
time_step = 10.0 # time step length
N_communities = 1 # how many communites there are in the population
community_probs = [1.0] # probability of being in each community
community_contact_rate = [1.0]
density_dependent_fecundity = 0.0007
average_worm_lifespan = 5.7
max_age = 100
initial_worms = 0
initial_miracidia_days = 3
initial_miracidia = 50000*N/1000
init_env_cercariae = 50000*N/1000
worm_stages = 1

female_factor = 1
male_factor = 1
age_contact_rates = [0.032,0.61, 1,0.06]
ages_for_contacts = [4,9,15,100]
contact_rate_by_age_array = fill(0.09, trunc(Int,(max_age+1)))
mda_adherence = 0.9
mda_access = 0.9
female_factor = 1
male_factor = 1
birth_rate = 28*time_step/(1000*365)
human_cercariae_prop = 1
predis_aggregation = 0.24
cercariae_survival = 1/2
miracidia_survival = 1/2
miracidia_maturity = 24
death_prob_by_age = [0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                     0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1]

ages_for_death = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                   65, 70, 75, 80, 85, 90, 95, 100, 110];
r = 0.03
kato_katz_par = 0.87
heavy_burden_threshold = 16
vaccine_effectiveness  = 0.86
drug_effectiveness = 0.86
spec_ages = [7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
              2636, 1970, 1468, 1166, 943, 718, 455, 244]
ages_per_index = 5
record_frequency = 1/24;

scenario = "moderate adult";
filename = "eq_runs.jld";
use_kato_katz = 1
# main parameters that we change
predis_aggregation = 0.24
contact_rate = 0.16
max_fecundity = 0.14
