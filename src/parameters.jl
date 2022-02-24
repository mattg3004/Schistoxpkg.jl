using JLD
using Distributions
using Random

N = 1  #population size
time_step = 1.0 # time step length
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
egg_production_distribution = "NegBin"
female_factor = 1
male_factor = 1
age_contact_rates = [0.1094983,
0.4363393,
0.4291753,
0.02498713]
age_contact_rates = age_contact_rates/sum(age_contact_rates)
c1 = age_contact_rates[1]
c2 = age_contact_rates[2]
c3 = age_contact_rates[3]
c4 = age_contact_rates[4]
ages_for_contacts = [4,9,15,100]
contact_rate_by_age_array = fill(0.00, trunc(Int,(max_age+1)))
mda_adherence = 1
mda_access = 1

last_uptake = 0
egg_multiplier = 100
sd_decrease = 1
female_factor = 1
male_factor = 1
birth_rate = 28*time_step/(1000*365)
human_cercariae_prop = 1
predis_aggregation = 0.24
cercariae_survival = 0.1
miracidia_survival = 0.1
miracidia_maturity = 24
death_prob_by_age = [0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                     0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1]

ages_for_death = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                   65, 70, 75, 80, 85, 90, 95, 100, 110];
r = 0.03
kato_katz_par = 0.87

vaccine_effectiveness  = 0.86
drug_effectiveness = 0.86
spec_ages = [7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
              2636, 1970, 1468, 1166, 943, 718, 455, 244]
ages_per_index = 5
record_frequency = 1/24;

scenario = "moderate adult";
filename = "eq_runs.jld";
use_kato_katz = 0
# main parameters that we change

predis_aggregation = 0.08489151
max_fecundity = 100
max_fec_contact_rate_product = 1/4
contact_rate = max_fec_contact_rate_product/max_fecundity

M0 = 1000
rate_acquired_immunity = 0
human_larvae_maturity_time = 40
egg_sample_size = 1/100
heavy_burden_threshold = 400


pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
        worm_stages, contact_rate, max_fec_contact_rate_product, max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
        death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index, record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
        rate_acquired_immunity, M0, human_larvae_maturity_time, egg_sample_size, egg_production_distribution)
