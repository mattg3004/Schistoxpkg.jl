using Schistoxpkg
using Test
using Distributions
using Random


env = Environment()
humans = env.humans

N = 1000    #population size
time_step = 10.0
N_communities = 1
community_probs = [1.0]
community_contact_rate = [1.0]
density_dependent_fecundity = 0.0007
average_worm_lifespan = 5.7
max_age = 100
initial_worms = 0
initial_miracidia = 100000*N/1000
initial_miracidia_days = 3
init_env_cercariae = 100000*N/1000
worm_stages = 1
contact_rate = 0.1
max_fecundity = 0.34
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
vaccine_effectiveness  = 0.86
drug_effectiveness = 0.86
spec_ages = [7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
              2636, 1970, 1468, 1166, 943, 718, 455, 244]
ages_per_index = 5
record_frequency = 1/24;
predis_aggregation = 0.8
scenario = "moderate adult";
use_kato_katz = 0
pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
        worm_stages, contact_rate, max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
        death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index, record_frequency, use_kato_katz)
pars = make_age_contact_rate_array(pars, scenario, [],[]);

humans, miracidia, cercariae = create_population_specified_ages(pars)

@testset "create_population" begin
    @test length(humans) == N
end



humans = generate_ages_and_deaths(20000, humans, pars)
humans = update_contact_rate(humans,  pars)
@testset "update_contact_rate1" begin
    @test humans[1].age_contact_rate == pars.contact_rate_by_age_array[end]
end
humans[1].age = 0
humans = update_contact_rate(humans,  pars)

@testset "update_contact_rate2" begin
    @test humans[1].age_contact_rate == pars.contact_rate_by_age_array[1]
end

@testset "miracidia_death" begin
    @test isapprox(miracidia_death!(miracidia, pars), [pars.initial_miracidia, pars.initial_miracidia , round(pars.initial_miracidia  * pars.miracidia_survival)])
end


@testset "cercariae_death" begin
    @test cercariae_death!(cercariae, pars)== pars.init_env_cercariae  * pars.cercariae_survival

end



@testset "worm_pairs" begin
    @test isapprox(calculate_worm_pairs([[1,2],[2,2]],[[5,1],[2,0]]), [3,2])
end

humans[1].female_worms = [1,2]
humans[1].male_worms = [5,1]

female_worms = (p->p.female_worms).(humans)
male_worms = (p->p.male_worms).(humans)

@testset "worm_pairs" begin
    @test calculate_worm_pairs(female_worms, male_worms)[1] ==3
end


@testset "make_age_contact_rate_array(max_age,scenario)" begin
    @test pars.contact_rate_by_age_array[1] == 0.032
end


@testset "make_age_contact_rate_array(max_age,scenario)" begin
    @test pars.contact_rate_by_age_array[100] == 0.06
end
