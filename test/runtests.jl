using Schistoxpkg
using Test
using Distributions
using Random


# Test that at the first element of the death array created
# matches the equation we think it should match
@testset "make_death_rate_array" begin
    @test make_death_rate_array([6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48, 0.53, 0.65,
    0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
    21.83, 29.98, 36.98], 1)[1] == 1 .- exp(-1 * 6.56/(1000*365))
end

# Test that the outputs of death rate array are sensible
death_array1 = make_death_rate_array([6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 0.88, 1.06,
1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 21.83, 29.98, 36.98], 1)

# We put in 21 parameters, so should be a length 21 array
@testset "make_death_rate_array_size" begin
    @test size(death_array1)[1] == 21
end

# They should all be positive.
@testset "make_death_rate_array_positive" begin
    @test all(death_array1 .> 0)
end




death_rate = find_death_rate(0, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])
@testset "find_death_rate" begin
    @test death_rate == 1
end

death_rate_2 = find_death_rate(6, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
@testset "find_death_rate2" begin
    @test death_rate_2 == 3
end

death_rate_3 = find_death_rate(1000, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
@testset "find_death_rate3" begin
    @test death_rate_3 == 18
end

# If I understand right, death rate should increase with age.
@testset "find_death_rate3" begin
    @test death_rate < death_rate_2
end

@testset "update_contact_rate" begin
    @test isapprox(update_contact_rate([5,10,3], [0,0,0], [1,2,3,4,5,6,7,8,9,10,11,12]), [6,11,4])
end



@testset "miracidia_death" begin
    @test miracidia_death([1000,1421,43523,55],1/2) == [1000,1421,43523,28]
end


@testset "miracidia_death" begin
    @test miracidia_death([1000,1421,43523,55],1/3) == [1000,1421,43523,18]
end



@testset "cercariae_death" begin
    @test cercariae_death(1000,1/2,1) == 500
end


@testset "cercariae_death" begin
    @test cercariae_death(1000,1/3,1) == 333
end


@testset "worm_pairs" begin
    @test isapprox(calculate_worm_pairs([[1,3],[17,0],[400,12312]],[[100,3000],[1,6],[1,1]]),[4,7,2] )
end

@testset "total_worms" begin
    @test isapprox(calculate_total_worms([[1,3],[17,0],[400,1]],[[10,30],[1,6],[1,1]])[1], [4,17,401])
end


@testset "total_worms" begin
    @test isapprox(calculate_total_worms([[1,3],[17,0],[400,1]],[[10,30],[1,6],[1,1]])[2], [40,7,2])
end



@testset "update_death_rate" begin
    @test update_death_rate([1,4,55,8], [1,1,1,1], [1,2,3,4,5,6,7,8,9,10,11,12,13]) == [2,2,12,3]
end



@testset "make_age_contact_rate_array(max_age,scenario)" begin
    @test make_age_contact_rate_array(100, "high adult", [], [])[1] == 0.01
end

@testset "make_age_contact_rate_array(max_age, scenario)" begin
    @test make_age_contact_rate_array(100, "high adult", [], [])[100] == 0.12
end

@testset "make_age_contact_rate_array(max_age, scenario)" begin
    @test make_age_contact_rate_array(100, "low adult", [], [])[100] == 0.02
end

@testset "make_age_contact_rate_array(max_age, scenario)" begin
    @test make_age_contact_rate_array(100, "moderate adult", [], [])[100] == 0.06
end
time_step = 10
human_cercariae = [[4],[5],[6],[7],[8]]
env_miracidia = [10, 10, 10, 10, 1]
env_cercariae = 100
contact_rate = 0.01
predisposition  = [0.1, 0.21, 1.4, 3.1, 100]
age_contact_rate = [0.0, 0.06, 0.12, 0.02, 0.06]
vac_status = [0,0,0,0,0]
vaccine_effectiveness = 0.95
community = 1,2,3,4,1
community_contact_rate = 1,1,1,0,1
female_worms= [[1,1], [2,1], [2,4],[3,1], [2,4]]
male_worms= [[1,1], [2,1], [2,5],[2,1], [2,4]]

@testset "cercariae_uptake(max_age, scenario)" begin
    @test isapprox(cercariae_uptake(copy( env_miracidia), copy(env_cercariae), copy(time_step), copy(contact_rate),
    community, community_contact_rate, copy(female_worms), copy(male_worms),
        copy(predisposition), copy(age_contact_rate), copy(vac_status), copy(vaccine_effectiveness), 1, 24)[2],  [10, 10, 10, 1])
end

time_step = 10
human_cercariae = [[4],[5],[6],[7],[8]]
env_miracidia = [10, 10, 10, 10, 1]
env_cercariae = 100
contact_rate = 0.01
predisposition  = [0.1, 0.21, 1.4, 3.1, 100]
age_contact_rate = [0.0, 0.06, 0.12, 0.02, 0.06]
vac_status = [0,0,0,0,0]
vaccine_effectiveness = 0.95

@testset "cercariae_uptake" begin
    @test isapprox(cercariae_uptake(copy( env_miracidia), copy(env_cercariae), copy(time_step), copy(contact_rate),
    community, community_contact_rate, copy(female_worms), copy(male_worms),
        copy(predisposition), copy(age_contact_rate), copy(vac_status), copy(vaccine_effectiveness),1,24)[3][1][2],  1)
end


time_step = 10
human_cercariae = [[4],[5],[6],[7],[8]]
env_miracidia = [10, 10, 10, 10, 1]
env_cercariae = 100
contact_rate = 0.01
predisposition  = [0.1, 0.21, 1.4, 3.1, 100]
age_contact_rate = [0.0, 0.06, 0.12, 0.02, 0.06]
vac_status = [0,0,0,0,0]
vaccine_effectiveness = 0.95

@testset "cercariae_uptake" begin
    @test cercariae_uptake(copy( env_miracidia), copy(env_cercariae), copy(time_step), copy(contact_rate),
    community, community_contact_rate, copy(female_worms), copy(male_worms),
        copy(predisposition), copy(age_contact_rate), copy(vac_status), copy(vaccine_effectiveness),1,24)[3][5][2] >0
end


time_step = 10
human_cercariae = [[4],[5],[6],[7],[8]]
env_miracidia = [10, 10, 10, 10, 1]
env_cercariae = 100
contact_rate = 0.01
predisposition  = [0.1, 0.21, 1.4, 3.1, 100]
age_contact_rate = [0.0, 0.06, 0.12, 0.02, 0.06]
vac_status = [0,0,0,0,1]
vaccine_effectiveness = 1

@testset "cercariae_uptake" begin
    @test cercariae_uptake(copy( env_miracidia), copy(env_cercariae), copy(time_step), copy(contact_rate),
    community, community_contact_rate, copy(female_worms), copy(male_worms),
        copy(predisposition), copy(age_contact_rate), copy(vac_status), copy(vaccine_effectiveness),1,1)[3][5][2] == 4
end

female_worms = [[10000,200000],[200000,0]]
male_worms = [[10000,200000],[200000,0]]
worm_stages = 2
average_worm_lifespan = 5.7

@testset "worm_maturity" begin
    @test worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)[1][1][2] < 200000
end

@testset "worm_maturity" begin
    @test worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)[1][2][2] > 00
end


@testset "worm_maturity" begin
    @test worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)[2][1][1] < 10000
end

@testset "worm_maturity" begin
    @test worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)[2][2][1] < 200000
end



@testset "miracidia_production" begin
    @test miracidia_production([1,2,3],[2,4,1], 10, [1,1,1],[1,1,2], [1,2,3])[1] == 2
end


@testset "miracidia_production" begin
    @test miracidia_production([1,2,3],[2,4,1], 10, [1,1,1],[1,1,1.2], [1,2,3])[end] == 6
end


@testset "miracidia_production" begin
    @test miracidia_production([1,2,3],[2,4,1], 10, [1,1,1],[0,1,0], [1,2,3])[end] == 2
end


@testset "miracidia_production" begin
    @test miracidia_production([1,2,3],[2,4,1], 10, [1,1,1],[0.5,0.5,0.5], [1,2,3])[end] == 6
end






vaccine_info = []
push!(vaccine_info, vaccine_information(0.75, 4, 16, [0,1], 3, 13))
push!(vaccine_info, vaccine_information(0.4, 17, 110, [0,1], 3, 17))

@testset "update_vaccine" begin
    @test update_vaccine(vaccine_info, 1)[1] == 0.4
end

@testset "update_vaccine" begin
    @test update_vaccine(vaccine_info, 1)[2] == 17
end

@testset "update_vaccine" begin
    @test update_vaccine(vaccine_info, 1)[3] ==110
end

@testset "update_vaccine" begin
    @test update_vaccine(vaccine_info, 1)[4] ==17
end

@testset "update_vaccine" begin
    @test update_vaccine(vaccine_info, 1)[5] ==[0,1]
end


@testset "update_vaccine" begin
    @test update_vaccine(vaccine_info, 100)[4] == Inf
end


@testset "death_of_human" begin
    @test death_of_human([2,4], [0,6], [1,0], [0.0002,0.00005], [[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1,1], 1,
                                [1,1], [1,1])[1] == [4]

end

@testset "death_of_human" begin
    @test death_of_human([2,4], [0,6], [1,0], [0.0002,0.00005], [[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1,1], 1,
                                [1,1], [1,1])[4] == [0.00005]

end


@testset "death_of_human" begin
    @test death_of_human([120,120], [0,0], [0.4,0.6], [[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1,1], [0.00001,0.0001], 365,
                                [1,1], [1,1])[1] == []

end
age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                           21.83, 29.98, 36.98]

contact_rates_by_age = make_age_contact_rate_array(100,"high adult", [], [])
death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, 1)

#
# ages, death_ages, gender, predisposition, human_cercariae, eggs, vac_status,
#                         treated, female_worms, male_worms,vaccinated, age_contact_rate,
#                         female_factor, male_factor, contact_rates_by_age,
#                         worm_stages, predis_aggregation, predis_weight,
#                         adherence, death_prob_by_age, ages_for_deaths,
#                         mda_adherence, access, mda_access

death_prob_by_age = [0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                          0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1]

ages_for_deaths = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                        65, 70, 75, 80, 85, 90, 95, 100, 110]



    @testset "birth_of_human" begin
        @test birth_of_human([2,4], [0.44,0], [0,0], [1.1,1], [2,4],[[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1.1,1], 1, 1, [0.0002,0.00005], 2,1, 1, [1,1],
                            death_prob_by_age, ages_for_deaths,[1,2,3,4], 0.8,[1,1], 0.9)[1] == [2,4,0]
    end
    new_pop = birth_of_human([2,4], [0.44,0], [0,0], [1.1,1],[2,4], [[2,3,4],[6,3,4]], [15,7],
                            [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                            [0,0], [1.1,1], 1, 1, [0.0002,0.00005], 2,1, 1, [1,1],
                        death_prob_by_age, ages_for_deaths,[1,2,3,4,5], 0.8,[1,1], 0.9)
    @testset "birth_of_human" begin
        @test new_pop[1]==[2,4,0]
    end
    @testset "birth_of_human" begin
        @test new_pop[6]==[[2,3,4],[6,3,4],[]]
    end

    @testset "birth_of_human" begin
        @test new_pop[7]==[15,7,0]
    end



@testset "administer_drug" begin
    @test administer_drug([[3,2],[3],[4,10]], [1,2], 1, [1,1,1]) == [[0,0],[0],[4,10]]
end




@testset "administer_drug" begin
    @test administer_drug([[3,2],[3],[4,10]], [1,2], 0, [1,1,1]) == [[3,2],[3],[4,10]]
end

Random.seed!(2525)
big_drug = administer_drug([[3, 2],[1e5],[4,10]], [2], 0.5, [1,1,1])
# I ran it with this seed so know the exact answer, so consider this a regression test.
@testset "administer_drug" begin
    @test big_drug == [[3,2],[50029],[4,10]]
end


@testset "create_mda" begin
     mda_info = create_mda(0, .75, 1, 1, 5, 2, [0,1], [0,1], [0,1], .92)
    @test mda_info[1].coverage == 0
end


@testset "create_mda" begin
     mda_info = create_mda(0, .75, 1, 1, 5, 1, [0,1], [0,1], [0,1], .92)
    @test mda_info[end].time == 5
end

mda_info = create_mda(0, .75, 1, 1, 5, 2, [0,1], [0,1], [0,1], .92)

@testset "update_mda" begin
    @test update_mda(mda_info, 3)[1] == 0
end

@testset "update_mda" begin
    @test update_mda(mda_info, 3)[5] == 3
end

@testset "update_mda" begin
    @test update_mda(mda_info, 30)[5] == Inf
end

  # But also with this big a sample size I think its fairly reasonable to use the 0.0001 and 0.9999 quantiles.
# In R qbinom(0.5, 1e5, c(0.0001, 0.9999))
# So this isn't relying on the seed, and should nearly always be true.
# but if we're repeatedly running it, we still want to use a seed.
@testset "administer_drug" begin
    @test (big_drug[2][1] > 10.0) & (big_drug[2][1] < 99990.0)
end


@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],0,0,[1,1,1,1,1], [0,0,0,0,0],[0,0,0,0,0])[4],
[[1], [2], [3], [4], [5]])
end



@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],0,0,[1,1,1,1,1], [1,1,1,1,1],[1,1,1,1,1])[1],
[[1], [0], [0], [0], [5]])
end

@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],0,0,[1,1,1,1,1], [1,1,1,1,1],[1,1,1,1,1])[2],
[[1], [0], [0], [0], [5]])
end

@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
0,0,[1,1,1,1,1], [1,1,1,1,1],[1,1,1,1,1])[3],
[[1], [0], [0], [0], [5]])
end


#=
Run some tests on the population making function
Need to do some setup first though.
=#
 mda_access = 1
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
#const contact_rate = 0.000005
age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                           21.83, 29.98, 36.98]
predis_aggregation = 0.24
mda_adherence = 0.8
scenario = "high adult"
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [], [])
death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)
predis_weight = 1

N_communities= 4
community_probs = [1,2,1,3]

pop = create_population(N, max_age, N_communities, community_probs, initial_worms, contact_rates_by_age,
    death_rate_per_time_step,worm_stages, female_factor, male_factor,
    initial_miracidia, initial_miracidia_days, predis_aggregation, predis_weight, time_step,
    mda_adherence, mda_access)

# These should all get give the initial value
@testset "miracidia" begin
    @test all(pop[14] .== initial_miracidia)
end

# Worms should be 0:ininity
# Get first column of worms
mworm1 = [pop[9][i][1] for i=1:length(pop[9])]
fworm1 = [pop[10][i][1] for i=1:length(pop[10])]
@testset "wormsm" begin
    @test all(mworm1 .>= 0)
end
@testset "wormsf" begin
    @test all(fworm1 .>= 0)
end

# All vectors except the last should be length N?
lens = [length(pop[i]) for i=1:(length(pop) - 3)]
push!(lens, length(pop[end]))
push!(lens, length(pop[end-1]))
@testset "N" begin
    @test all(lens .== N)
end





female_worms= [[1,1], [2,1], [2,4]]
male_worms= [[1,1], [2,1], [2,5]]
vaccine_effectiveness = 1
vaccine_coverage = 1
min_age_vaccine = 5
max_age_vaccine = 16
vaccine_gender = [0,1]
ages = [4,10,9]
death_ages = [5,11,66]
human_cercariae = [[1,2],[9,2,4], [ 0,5]]
eggs = [1,2,3]
treated = [1,1,1]
vaccine_duration = 10
vac_status = [0,0,0]
vaccine_round = 1
gender = [1,0,1]
adherence = [1,1,1]
access = [1,1,1]


@testset "vaccinate" begin
@test isapprox(vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
        vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
        treated, vaccine_duration, vac_status, vaccine_round, gender, access)[1],
        [[1,1],[0,0],[0,0]])
end

female_worms= [[1,1], [2,1], [2,4]]
male_worms= [[1,1], [2,1], [2,5]]
access = [1,1,1]



@testset "vaccinate" begin
@test isapprox(vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
            vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
            treated, vaccine_duration, vac_status, vaccine_round, gender,access)[1],
            [[1,1],[0,0],[0,0]])
end

female_worms= [[1,1], [2,1], [2,4]]
male_worms= [[1,1], [2,1], [2,5]]
access = [1,0,0]

@testset "vaccinate" begin
@test isapprox(vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
            vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
            treated, vaccine_duration, vac_status, vaccine_round, gender,access)[1],
            [[1,1],[2,1],[2,4]])
end

time_step = 10
num_time_steps = 1
community = 1,3,1
community_contact_rate = 1,0.5, 1
community_probs = 1,2,1

@testset "update_env" begin
    @test update_env(num_time_steps, [1,3,4], [2,5,8],community, community_contact_rate, community_probs, human_cercariae, female_worms, male_worms,
    time_step, 5.7,
    eggs, 0.34, 0.03, 2,
    vac_status, gender, 0.24,predis_weight,
    [0.2, 0.8,1], treated, vaccine_effectiveness,
    0.0005, death_prob_by_age, ages_for_deaths,
    [0,0,0], [0.02,0.04,0.03], env_miracidia,
    env_cercariae, 0.0005, 1, 1,
    1, 1, contact_rates_by_age,
    28*time_step/(1000*365), [], [], [1,1,1], 1,
    [1,1,1], 1,
    1/24,1, 24)[1] ==  [1+(num_time_steps*time_step/365),3+(num_time_steps*time_step/365),4+(num_time_steps*time_step/365),24]
end
