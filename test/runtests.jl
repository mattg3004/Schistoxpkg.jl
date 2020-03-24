using Schistoxpkg

using Test
using Distributions
using Random


# Test that at the first element of the death array created
# matches the equation we think it should match
@testset "make_death_rate_array" begin
    @test make_death_rate_array([6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
    0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
    21.83, 29.98, 36.98], 1)[1] == 1 .- exp(-1 * 6.56/(1000*365))
end

# Test that the outputs of death rate array are sensible
death_array1 = make_death_rate_array([6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 21.83, 29.98, 36.98], 1)

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





@testset "update_death_rate" begin
    @test isapprox(update_death_rate([0,3,6,50, 1000], [0,0,0,0,0], [1,2,3,4,5,6,7,8,9,10,11,12]), [1,2,3,11,12])
end



@testset "cercariae_death" begin
    @test cercariae_death(1000,1,1) == 500
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


@testset "make_age_contact_rate_array(max_age)" begin
    @test make_age_contact_rate_array(100)[1] == 0.032
end




@testset "death_of_human" begin
    @test death_of_human([2,4], [0,0], [0.4,0.6], [[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1,1], [0.0002,0.00005], 1,
                                [1,1])[1] == [2,4]

end

@testset "death_of_human" begin
    @test death_of_human([2,4], [0,0], [0.4,0.6], [[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1,1], [2,0.00005], 1,
                                [1,1])[1] == [4]

end


@testset "death_of_human" begin
    @test death_of_human([120,120], [0,0], [0.4,0.6], [[2,3,4],[6,3,4]], [15,7],
                                [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                                [0,0], [1,1], [0.00001,0.0001], 365,
                                [1,1])[1] == []

end
age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                           21.83, 29.98, 36.98]

contact_rates_by_age = make_age_contact_rate_array(100)
death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, 1)

@testset "birth_of_human" begin
    @test birth_of_human([2,4], [0,0], [0.4,0.6], [[2,3,4],[6,3,4]], [15,7],
                            [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                            [0,0], [1.1,1], [0.0002,0.00005], 1,1, contact_rates_by_age,
                        death_rate_per_time_step,2, 0.24, [1,1],1)[1] == [2,4,0]
end
new_pop = birth_of_human([2,4], [0,0], [0.4,0.6], [[2,3,4],[6,3,4]], [15,7],
                        [0,0], [0,0], [[9,2],[5,3]], [[0,3],[1,12]],
                        [0,0], [1.1,1], [0.0002,0.00005], 1,1, contact_rates_by_age,
                    death_rate_per_time_step,2, 0.24, [1,1],1)
@testset "birth_of_human" begin
    @test new_pop[1]==[2,4,0]
end
@testset "birth_of_human" begin
    @test new_pop[4]==[[2,3,4],[6,3,4],[]]
end

@testset "birth_of_human" begin
    @test new_pop[5]==[15,7,0]
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
     mda_info = create_mda(0, .75, 1, 1,5, 2, [0,1], [0,1], [0,1], .92)
    @test mda_info[1].coverage == 0
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
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],0,0,[1,1,1,1,1], [0,0,0,0,0])[4],
[[1], [2], [3], [4], [5]])
end



@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],0,0,[1,1,1,1,1], [1,1,1,1,1])[1],
[[1], [0], [0], [0], [5]])
end

@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],0,0,[1,1,1,1,1], [1,1,1,1,1])[2],
[[1], [0], [0], [0], [5]])
end

@testset "mda" begin
    @test isapprox(mda(1, 3, 8, 1,[0,1],
[2,5,4,7,9], [[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
[[1],[2],[3],[4],[5]], [[1],[2],[3],[4],[5]],
0,0,[1,1,1,1,1], [1,1,1,1,1])[3],
[[1], [0], [0], [0], [5]])
end


#=
Run some tests on the population making function
Need to do some setup first though.
=#

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

contact_rates_by_age = make_age_contact_rate_array(max_age)
death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)

pop = create_population(N, max_age, initial_worms, contact_rates_by_age,
    death_rate_per_time_step,worm_stages, female_factor, male_factor,
    initial_miracidia, initial_miracidia_days, predis_aggregation, time_step,
    mda_adherence)

# These should all get give the initial value
@testset "miracidia" begin
    @test all(pop[13] .== initial_miracidia)
end

# Worms should be 0:ininity
# Get first column of worms
mworm1 = [pop[8][i][1] for i=1:length(pop[8])]
fworm1 = [pop[9][i][1] for i=1:length(pop[9])]
@testset "wormsm" begin
    @test all(mworm1 .>= 0)
end
@testset "wormsf" begin
    @test all(fworm1 .>= 0)
end

# All vectors except the last should be length N?
lens = [length(pop[i]) for i=1:(length(pop) - 2)]
push!(lens, length(pop[end]))
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
human_cercariae = [[1,2],[9,2,4], [ 0,5]]
eggs = [1,2,3]
treated = [1,1,1]
vaccine_duration = 10
vac_status = [0,0,0]
vaccine_round = 1
gender = [1,0,1]
adherence = [1,1,1]


@testset "vaccinate" begin
@test isapprox(vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
        vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
        treated, vaccine_duration, vac_status, vaccine_round, gender, adherence)[1],
        [[1,1],[0,0],[0,0]])
end

female_worms= [[1,1], [2,1], [2,4]]
male_worms= [[1,1], [2,1], [2,5]]
adherence = [1,0,0]

@testset "vaccinate" begin
@test isapprox(vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
            vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
            treated, vaccine_duration, vac_status, vaccine_round, gender, adherence)[1],
            [[1,1],[2,1],[2,4]])
end
