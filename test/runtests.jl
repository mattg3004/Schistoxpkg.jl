using Schistoxpkg
using Test
using Distributions
using Random
using JLD



egg_production_distribution = "NegBin"
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
max_fec_contact_rate_product =  0.034
contact_rate = max_fec_contact_rate_product / max_fecundity
female_factor = 1
male_factor = 1
age_contact_rates = [0.032,0.61, 1,0.06]
ages_for_contacts = [4,9,15,100]
contact_rate_by_age_array = fill(0.09, trunc(Int,(max_age+1)))
mda_adherence = 1
mda_access = 1
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
kato_katz_par = 0.87
heavy_burden_threshold = 16
rate_acquired_immunity = 0
M0 = 15
egg_sample_size  = 1
filename = "ee.jld"
human_larvae_maturity_time = 40



@testset "test_pars" begin
    pars = Parameters()
    @test pars.N == 1000
end


pars = Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
        worm_stages, contact_rate, max_fec_contact_rate_product,  max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
        death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index, record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold, rate_acquired_immunity, M0,
        human_larvae_maturity_time, egg_sample_size, egg_production_distribution)
pars = make_age_contact_rate_array(pars, scenario, [],[]);

humans, miracidia, cercariae = create_population_specified_ages(pars)
@testset "test_human" begin
    h = Human()
    @test typeof(h) == Human
end




@testset "create_contacts_settings" begin
    conts = create_contact_settings("low adult")
    @test conts[4] == 0.02
end


@testset "create_contacts_settings" begin
    conts = create_contact_settings("high adult")
    @test conts[4] == 0.12
end


@testset "create_contacts_settings" begin
    conts = create_contact_settings("moderate adult")
    @test conts[4] == 0.06
end

pars.max_age = 5


@test_throws ErrorException("max_age must be greater than 60") make_age_contact_rate_array(pars, "high adult", [], [])
pars.max_age = 100
pars = make_age_contact_rate_array(pars, "high adult", [], [])
@testset "contact_array" begin

    @test pars.contact_rate_by_age_array[1] == 0.01
end


humans, miracidia, cercariae = create_population_specified_ages(pars)

@testset "create_population" begin
    @test length(humans) == N
end



humans, miracidia, cercariae = create_population(pars)

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


updated_cercariae = cercariae
time_step_specific_cerc_death = pars.cercariae_survival
c = cercariae
for i in 1:pars.time_step
    c = (c + miracidia[1]/pars.time_step) * pars.cercariae_survival
end
c1 = cercariae
c1 = (c1 +  miracidia[1])
updated_cercariae = trunc(Int, round(cercariae * time_step_specific_cerc_death, digits= 0))


@testset "cercariae_death" begin
    @test cercariae_death!(cercariae, miracidia, pars)== updated_cercariae

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
    @test pars.contact_rate_by_age_array[100] == 0.12
end


humans, miracidia, cercariae = create_population(pars)
old_cercariae = cercariae
old_length_miracidia = length(miracidia)



humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)
@testset "cerc_uptake" begin
    @test length(miracidia) <  old_length_miracidia
end

human_larvae = (p->p.larvae).(humans)



@testset "cerc_uptake" begin
    @test length(miracidia) <  old_length_miracidia
end

@testset "human_larvae" begin
    @test sum(sum.(human_larvae)) == 0
end

cercariae = 10000
humans, cercariae, miracidia  = cercariae_uptake_with_human_larvae!(humans, cercariae, miracidia, pars)


human_larvae = (p->p.larvae).(humans)

@testset "human_larvae" begin
    @test sum(sum.(human_larvae)) > 0
end

human_larvae = (p->p.larvae).(humans)
x = sum(sum.(human_larvae))
humans = human_larvae_maturity(humans, pars)
human_larvae = (p->p.larvae).(humans)

@testset "human_larvae_maturity" begin
    @test sum(sum.(human_larvae)) == x
end

pars = make_age_contact_rate_array(pars, scenario, [4,9,15,pars.max_age], [10.1,2.1,4.2,6.2])
@testset "new_contact_rate" begin
    @test pars.contact_rate_by_age_array[1] == 10.1
end
start_worms = 10000
for h in humans
    h.female_worms = [start_worms]
    h.male_worms = [start_worms]
end
humans = worm_maturity!(humans, pars)
@testset "worm_maturity" begin
    @test humans[1].female_worms[1] < start_worms
end

#
# @testset "calc_worm_pairs" begin
#     @test isapprox(calculate_worm_pairs([[1,1],[2,3]], [[2,3],[5,3]]) , [2,5])
# end

pars.mda_adherence = 1
pars.mda_access = 1
pars.initial_worms = 50
Random.seed!(2525)
humans, miracidia, cercariae = create_population(pars)


humans = egg_production!(humans, pars)
# @testset "egg_prod" begin
#     @test mean((p->p.eggs).(humans)) == 5.402
# end
#
# @testset "mira_prod" begin
#     @test miracidia_production!(humans) == 81
# end
eggs = sum((p->p.eggs).(humans))
humans = death_of_human(humans)
@testset "death" begin
    length(humans) < pars.N
end


humans = administer_drug(humans, 1:length(humans), 1)
@testset "administer_drug" begin
    @test sum((p->p.eggs).(humans)) == eggs
end


@testset "administer_drug" begin
    @test sum(sum.((p->p.female_worms).(humans))) == 0
end


administer_vaccine(humans, 1:length(humans), 1, 3)

@testset "administer_vaccine" begin
    @test sum(sum.((p->p.vac_status).(humans))) == 3*length(humans)
end

Random.seed!(2525)
humans, miracidia, cercariae = create_population(pars)


humans = egg_production!(humans, pars)

mda(humans, 1, 0, 100, 1, [0,1])
@testset "mda" begin
    @test sum((p->p.eggs).(humans)) == 0
end


first_mda_time = 1
last_mda_time = 10
mda_info = create_mda(0, 0.75, 0, first_mda_time,
    last_mda_time, 1, [0,1], [0,1], [0,1], 1)
@testset "create_mda" begin
    @test length(mda_info) == (last_mda_time - first_mda_time)*3
end




mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender = update_mda(mda_info, 0)
@testset "update_mda" begin
    @test [min_age_mda, max_age_mda] == [0,4]
end


@testset "update_mda" begin
 @test update_mda(mda_info, 1)[1] == 0.75
end

@testset "vac_decay" begin
    @test vac_decay!(humans,pars)[1].vac_status == -pars.time_step
end

Random.seed!(33)
@testset "Kato_katz" begin
    @test kato_katz(10, Gamma(0.87, 1/0.87)) == 1
end

Random.seed!(33)
humans, miracidia, cercariae = create_population(pars)
humans = egg_production!(humans, pars)

# @testset "count_eggs" begin
#     @test count_eggs(humans) == 5597
# end
#
#
# record = get_prevalences!(humans, time, pars)
# @testset "get_prevs" begin
#     @test isapprox(record.sac_burden, [93.94, 57.58, 1.01])
# end
#
save(filename, "humans", humans,  "miracidia", miracidia, "cercariae", cercariae, "pars", pars)
#
# humans_loaded,  miracidia_loaded, cercariae_loaded, pars_loaded = load_population_from_file(filename)
#
# @testset "load_pop" begin
#     @test humans_loaded[1].age == humans[1].age
# end


Random.seed!(33)


# @testset "update_env_to_equ" begin
#     @test update_env_to_equilibrium(1, humans, miracidia, cercariae, pars)[2][3] == 41
# end
#
#
# Random.seed!(33)
# @testset "update_env_to_equ_inc" begin
#     @test update_env_to_equilibrium_increasing(1, humans, miracidia, cercariae, pars)[2][3] == 14
# end
#
vaccine_info = []
mda_info = create_mda(0, .75, 0, 0, 2, 1, [0,1], [0,1], [0,1], pars.drug_effectiveness)
Random.seed!(33)

#
# @testset "update_env_const_pop" begin
#     @test isapprox(update_env_constant_population(10, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2], [10160, 12146, 12620])
# end
#
#
#
#
# @testset "update_env_no_birth_death" begin
#     @test isapprox(update_env_no_births_deaths(10, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2], [19110, 19696, 20028])
# end
#
#
#
#
#
# @testset "update_env_human_larvae" begin
# @test isapprox(update_env_to_equilibrium_human_larvae(5, humans, miracidia, cercariae, pars)[2],[20672, 20682, 20594])
# end
#
#
# @testset "update_env_human_larvae_constant_pop" begin
#     @test isapprox(update_env_constant_population_human_larvae(5, humans,  miracidia,
#                     cercariae, pars, mda_info, vaccine_info)[2],
#                         [23944, 24494, 25000])
# end
#
#
#
#
# @testset "update_env_no_births_deaths_human_larvae" begin
# @test isapprox(update_env_no_births_deaths_human_larvae(5, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2],[27568, 28168, 28605])
# end
#
#
# @testset "update_env_no_births_deaths_increasing" begin
# @test isapprox(update_env_no_births_deaths_increasing(5, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2],[870, 908, 888])
# end
#
# @testset "update_env_constant_population_increasing" begin
# @test isapprox(update_env_constant_population_increasing(5, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2],[868, 888, 877])
# end

humans, miracidia, cercariae, record =
    update_env_constant_population(25, humans,  miracidia, cercariae, pars, mda_info, vaccine_info);


times = []
prev = []
sac_prev = []
high_burden = []
high_burden_sac =[]
adult_prev = []
high_adult_burden = []

times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden =
collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, 1)


times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden =
collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, 2)


#
# @testset "run_repeated_sims_no_population_change" begin
# @test run_repeated_sims_no_population_change(filename, 5, mda_info, vaccine_info, 1)[2][1][1] == 64.6
# end
#
# @testset "run_repeated_sims_no_population_change" begin
# @test run_repeated_sims_no_population_change_increasing(filename, 5, mda_info, vaccine_info, 1)[2][1][1] == 60.8
# end
#
# @testset "run_repeated_sims_no_population_change" begin
# @test run_repeated_sims_no_births_deaths(filename, 5, mda_info, vaccine_info, 1)[2][1][1] == 86.8
# end
#
# @testset "run_repeated_sims_no_births_deaths_human_larvae" begin
# @test run_repeated_sims_no_births_deaths_human_larvae(filename, 5, mda_info, vaccine_info, 1)[2][1][1] == 80.1
# end
#
#
# @testset "run_repeated_sims_no_births_deaths_increasing" begin
# @test run_repeated_sims_no_births_deaths_increasing(filename, 5, mda_info, vaccine_info, 1)[2][1][1] == 75.2
# end
