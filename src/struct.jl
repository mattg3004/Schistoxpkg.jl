using Plots
using Distributions

mutable struct Human
    age::Float32
    gender::Int64
    predisposition::Float32
    female_worms::Array{Int64}
    male_worms::Array{Int64}
    cercariae::Array{Int64}
    eggs::Int64
    vac_status::Int64
    treated::Array{Int64}
    vaccinated::Array{Int64}
    age_contact_rate::Float32
    death_age::Float32
    adherence::Int64
    access::Int64
end




mutable struct out
    population_burden::Array{Float32}
    sac_burden::Array{Float32}
    adult_burden::Array{Float32}
    pop_prev::Float32
    sac_prev::Float32
    adult_prev::Float32
    sac_pop::Int64
    adult_pop::Int64
    final_ages::Array{Float32}
    recorded_eggs::Array{Int64}
    time::Float32
end



# function to generate an age for death of an individual
function get_death_age(death_prob_by_age, ages_for_deaths)
        age = 0
        k = 1
        p = death_prob_by_age[k]
        goal_age = ages_for_deaths[k]
        x = rand()
        while x > p
            age += 1
            if age >= goal_age
                k += 1
                p = death_prob_by_age[k]
                goal_age = ages_for_deaths[k]
            end
            x = rand()
        end
        if k == 1
            min_age = 0
        else
            min_age = ages_for_deaths[k-1]
        end

        death_age = min_age + rand() * (goal_age - min_age)

    return death_age
end



# function to age population and generating death ages
function generate_ages_and_deaths(num_steps, humans, female_factor, male_factor, time_step,
                    death_prob_by_age, ages_for_deaths, worm_stages, predis_aggregation,
                    mda_adherence, mda_access)

     for i in 1:num_steps
        deaths = 0
        death_index = []
        for j in 1 : length(humans)
            if humans[j].age > humans[j].death_age
                deaths += 1
                push!(death_index, j)
            end
            humans[j].age += time_step/365
        end
        if deaths > 0
            for j in deaths:-1:1
                splice!(humans, death_index[j])
                humans = birth_of_human(humans, female_factor, male_factor,
                    death_prob_by_age, ages_for_deaths, worm_stages, predis_aggregation,
                    mda_adherence, mda_access)
                f_worms = fill(0, worm_stages)
                f_worms[1] = round(rand()*initial_worms)
                m_worms = fill(0, worm_stages)
                m_worms[1] = round(rand()*initial_worms)
                humans[end].female_worms = f_worms
                humans[end].male_worms = m_worms
            end
        end
    end
    return humans
end







# create the age specific contact settings given the scenario
function create_contact_settings(scenario)
    if scenario == "low adult"
        contact_settings = [0.01, 1.2, 1, 0.02]
    elseif scenario == "moderate adult"
        contact_settings = [0.032, 0.61, 1, 0.06]
    elseif scenario == "high adult"
        contact_settings = [0.01, 0.61, 1, 0.12]
    end
end


# function to get age dependent contact rate.
# the contact rates are taken from the
# "What is required in terms of mass drug administration to interrupt the transmission
#     of schistosome parasites in regions of endemic infection?" paper
# at some point we may change this to be an input from a file instead

function make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)
    if max_age < 60
        error("max_age must be greater than 60")
    else

        if length(input_ages) == 0
            contact_settings = create_contact_settings(scenario)
    # initialize an array with the same value for contact rate across all ages
            contact_rates_by_age = [fill(contact_settings[4], max_age+1)]
            contact_rates_by_age = contact_rates_by_age[1]

    # then edit the entries for different ages according to values
            for i in 1:5
                contact_rates_by_age[i] = contact_settings[1]
            end
            if scenario == "high adult"
                for i in 6:12
                    contact_rates_by_age[i] = contact_settings[2]
                end
                for i in 13:21
                    contact_rates_by_age[i] = contact_settings[3]
                end
            else
                for i in 6:10
                    contact_rates_by_age[i] = contact_settings[2]
                end
                for i in 11:16
                    contact_rates_by_age[i] = contact_settings[3]
                end
            end

        else
            contact_rates_by_age = [fill(input_contact_rates[end], max_age+1)]
            contact_rates_by_age = contact_rates_by_age[1]
            for i in 1 : length(input_contact_rates)
                if i == 1
                    for j in 1:(input_ages[i] + 1)
                        contact_rates_by_age[j] = input_contact_rates[i]
                    end
                else
                    for j in (input_ages[i-1]+2):(input_ages[i] + 1)
                        contact_rates_by_age[j] = input_contact_rates[i]
                    end
                end
            end

        end

        return contact_rates_by_age
    end
# return contact rates for age
end







# define a function which will create the initial population.
# This randomly chooses age, and gender.
# predisposition is taken to be gamma distributed.
# There is also a male and female adjustment to predisposition
# adjusting for gender specific behaviour

function create_population(N,
        max_age,
        initial_worms,
        age_contact_rates,
        predis_aggregation,
        male_factor,
        female_factor,
        initial_miracidia,
        initial_miracidia_days,
        worm_stages,
        death_prob_by_age,
        ages_for_deaths,
        mda_adherence,
        mda_access)

    #=  initialize and fill the environmental variable  =#
    env_miracidia = Int64[]

    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia )
    end

    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)
    humans = Human[]
    for i in 1:N
        if rand() > mda_adherence
            adherence = 0
        else
            adherence = 1
        end
        if rand() > mda_access
            access = 0
        else
            access = 1
        end
        f_worms = fill(0, worm_stages)
        f_worms[1] = round(rand()*initial_worms)
        m_worms = fill(0, worm_stages)
        m_worms[1] = round(rand()*initial_worms)
        push!(humans, Human(rand()*max_age, rand([0,1]), rand(gamma_pre)[1],
        f_worms, m_worms,
        [], 0, 0, [], [], 1, 0, adherence, access))
        age = trunc(Int, humans[end].age)
        humans[end].age_contact_rate = age_contact_rates[age+1]
        x = humans[end].gender
        humans[end].predisposition = humans[end].predisposition * (1-x) *female_factor + humans[end].predisposition * x *female_factor
        humans[end].death_age = get_death_age(death_prob_by_age, ages_for_deaths)
    end
    return humans, env_miracidia
end





#= function to generate a distribution for ages based on a specified demography =#
function generate_age_distribution(spec_ages, ages_per_index)

    number_per_age = []
    for i in 1:(length(spec_ages) * ages_per_index)
        index = trunc(Int, (i-1)/ages_per_index) + 1
        push!(number_per_age, spec_ages[index])
    end
    cumsum_spec_ages = cumsum(number_per_age)/sum(number_per_age)
    return cumsum_spec_ages
end



#= function to construct the set of ages, with size N =#
function specified_age_distribution(N, spec_ages, ages_per_index)

    cumsum_spec_ages = generate_age_distribution(spec_ages, ages_per_index)
    ages = Float32[]
    for i in 1:N
        x = rand()
        k = findall(cumsum_spec_ages .> x)[1]
        push!(ages, k-1)
    end
    return ages
end





function create_population_specified_ages(N, initial_worms, contact_rates_by_age,
        worm_stages, female_factor, male_factor,initial_miracidia,
        initial_miracidia_days, predis_aggregation, time_step,
        spec_ages, ages_per_index,
        mda_adherence, mda_access)


#=  initialize the Gamma distribution for predisposition selection  =#
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

#=  initialize and fill the environmental variable  =#
    env_miracidia = Int64[]

    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia )
    end

    predisposition = rand(gamma_pre,N)
    ages = specified_age_distribution(N, spec_ages, ages_per_index)
    humans = Human[]
    for i in 1:N

        if rand() > mda_adherence
            adherence = 0
        else
            adherence = 1
        end
        if rand() > mda_access
            access = 0
        else
            access = 1
        end
        f_worms = fill(0, worm_stages)
        f_worms[1] = round(rand()*initial_worms)
        m_worms = fill(0, worm_stages)
        m_worms[1] = round(rand()*initial_worms)
        push!(humans, Human(ages[i], rand([0,1]), predisposition[i],
        f_worms, m_worms,
        [], 0, 0, [], [], 1, 0, adherence, access))
        age = trunc(Int, humans[end].age)
        humans[end].age_contact_rate = age_contact_rates[age+1]
        x = humans[end].gender
        humans[end].predisposition = humans[end].predisposition * (1-x) *female_factor + humans[end].predisposition * x *female_factor
        humans[end].death_age = get_death_age(death_prob_by_age, ages_for_deaths)
    end

    return humans, env_miracidia
end





#= function to update the contact rate of individuals in the population =#

function update_contact_rate(d,  contact_rates_by_age)
    for i in 1:length(d)
        age = min(length(contact_rates_by_age) - 1, (trunc(Int, d[i].age)))
        @inbounds d[i].age_contact_rate =  contact_rates_by_age[age+1]
    end
    return d
end





# humans uptake larvae based on their predisposition, age_dependent contact rate
# and the number of larvae in the environment. number of larvae taken up is chosen
# from a Poisson distribution

function cercariae_uptake(d, env_miracidia, env_cercariae, time_step, contact_rate, vaccine_effectiveness)

#= we want the human population to randomly pick up larvae.
therefore we want to shuffle the population.
the following lines make a random permutation of indices for the population=#
    d = shuffle(d)
#= assign larvae which have been in the environment for 40 days to become infective.
then delete those larvae from the environmental larvae =#
    k = length(d)
    total_uptakes = 0
    if  length(env_miracidia) > (40 / time_step)
        env_cercariae += env_miracidia[1]
        splice!(env_miracidia, 1)
    end

#= loop over the population uptaking larvae=#
    for i in 1:length(d)

#= if there are still infective larvae in the environment,
we will uptake from choose from Poisson
distribution. otherwise, just uptake 0. =#
        if env_cercariae > 0
# calculate the rate of the poisson distribution
            @inbounds  pois_rate  = d[i].predisposition * contact_rate * d[i].age_contact_rate *
                    env_cercariae * time_step / k

# reduce the rate according to the effectiveness of the vaccine (if any is given)
            @inbounds  pois_rate = pois_rate * (1 - (d[i].vac_status > 0) * vaccine_effectiveness);

# choose from the Poisson distribution
            uptake = rand(Poisson(pois_rate))
    #println(uptake)
        else
            uptake = 0
        end
        total_uptakes += uptake
# push the number of uptaken larvae into the human larvae array
        push!(d[i].cercariae, uptake)

# # reduce the infective larvae by the number of larvae uptaken
        env_cercariae -= uptake

    end
#println(uptakes/k)
# return the infective, human and environmental larvae arrays
    return d, env_cercariae, env_miracidia, total_uptakes

end



# function to kill miracidia in the environment
function miracidia_death(env_miracidia, env_miracidia_survival_prop)
    #= as env_miracidia is an array, we need to use . syntax to apply
    the functions to each element in the array =#
    # return rand.(Binomial.(env_miracidia, 1 - env_miracidia_survival_prop))
    if env_miracidia_survival_prop <= 0
        error("env_miracidia_survival_prop must be bigger than 0")
    else
        env_miracidia[end] = trunc(Int, round(env_miracidia[end] * env_miracidia_survival_prop, digits = 0))
    end
    # for i in 1:length(env_miracidia)
    #     env_miracidia[i] = trunc(Int, round(env_miracidia[i]/1.5, digits = 0))
    # end
    return env_miracidia
end



# function to kill cercariae in the environment
function cercariae_death(env_cercariae, env_cercariae_survival_prop)
    # updated_cercariae = 0
    # for i in 1:time_step
    #     updated_cercariae += (env_cercariae/time_step) * (1/(2^i))
    # end
    if env_cercariae_survival_prop <= 0
        error("env_cercariae_survival_prop must be bigger than 0")
    else
        updated_cercariae = trunc(Int, round(env_cercariae * env_cercariae_survival_prop, digits= 0))
    end
    return updated_cercariae
end






# Worms die have a specified mean life span, and hence a rate of deaths per day .
# the number of deaths, or maturing from each stage is dependent
# on the number of worm stages and rate of aging through stages and dying.
# This p is multiplied by the time scale of the simulation, so if 2 days
# pass between consecutive time points, twice as many worms age and die


function worm_maturity(d, worm_stages,
    average_worm_lifespan, time_step)

    # probability of aging out of category/ dying
        p = time_step / ( worm_stages * 365 * average_worm_lifespan)



# loop over the worms
    for i in 1:length(d)


# kill appropriate number of worms in the final stage
        @inbounds n = d[i].female_worms[worm_stages]

        if n > 0
            dis = Binomial(n, 1-p)
            @inbounds d[i].female_worms[worm_stages] = rand(dis, 1)[1]
        end


        @inbounds n = d[i].male_worms[worm_stages]
        if n > 0
            dis = Binomial(n, 1-p)
            @inbounds d[i].male_worms[worm_stages] = rand(dis, 1)[1]
        end

#=
     for aging worms, we do this in reverse order, which ensures the
     correct order of aging is respected
=#

        for j in (worm_stages-1):-1:1
#=   choose the number of male and female worms to age from one stage to the next   =#
            @inbounds aging_females = rand(Binomial(d[i].female_worms[j], p), 1)[1]
            @inbounds aging_males = rand(Binomial(d[i].male_worms[j], p), 1)[1]

#=   add and subtract the number of worms from the appropriate categories   =#
            @inbounds d[i].female_worms[j+1] += aging_females
            @inbounds d[i].female_worms[j] -= aging_females
            @inbounds d[i].male_worms[j+1] += aging_males
            @inbounds d[i].male_worms[j] -= aging_males
        end
    end

#=  return the female and male worm arrays  =#
    return d

end







# function to calculate the number of male-female pairs of worms in each human
# this is just the minimum of total female worms and the total male worms


function calculate_worm_pairs(d)
    return min(sum(d.female_worms), sum(d.male_worms))
end




# function to calculate the number of male and female worms in each human


function calculate_total_worms(female_worms, male_worms)
    return sum(female_worms) + sum(male_worms)
end






#=
function to calculate the number of eggs produced
this is done by choosing from a negative binomial distribution for each worms,
where the mean and aggregation parameters are calculated as in the
"Refined stratified-worm-burden models that incorporate specific biological features
of human and snail hosts provide better estimates of Schistosoma diagnosis,
transmission, and control" paper
for julia the negative binomial describes the number of failures before
the given number of successes
in a collection of independent Bernoulli trials.
we need to specify a probability of success, and a given number of
successes, which are derived
from the mean and aggregation in the function below
=#

# inputs

# r - aggregation factor for NB distribution


function egg_production(d, max_fecundity, r, density_dependent_fecundity, time_step)


# loop over individuals
    for i in 1 : length(d)

#= if we have a positive number of worms, then make calculation,
otherwise the number of eggs is trivially 0 =#
        worm_pairs = calculate_worm_pairs(d[i])
        if worm_pairs > 0

# calculate the mean number of eggs we would expect
                 mean_eggs =  max_fecundity * worm_pairs *
                    exp(- density_dependent_fecundity *
                    (sum(d[i].female_worms)))

# calculate the number of successes
                NB_r = r * worm_pairs

# calculate the probability of a success
                p = NB_r/(NB_r + mean_eggs)

# choose from NB
                eggs_num = rand(NegativeBinomial(NB_r,p))[1]

            else
                eggs_num = 0
            end

# put this selected number of eggs into the eggs array
            @inbounds d[i].eggs = eggs_num
        end

# return the eggs array
    return d
end



# hatch the eggs in the humans into the environment

function miracidia_production(d, env_miracidia, age_contact_rates)
#= as we can step forward an arbitrary number of days at a time, we multiply the number of miracidia by the
    length of the forward step, assuming that each of the last given number of days were equivalent to each other
=#

#     max_contact_rate = maximum(age_contact_rates)
#     xx = age_contact_rate ./ max_contact_rate
#     released_eggs = xx .* eggs
#     push!(env_miracidia,  sum(released_eggs))
#     return env_miracidia

    max_contact_rate = maximum(age_contact_rates)
    new_miracidia = 0
    for i in 1:length(d)
         new_miracidia = new_miracidia + d[i].eggs * d[i].age_contact_rate / max_contact_rate
    end
    push!(env_miracidia, round(new_miracidia))
    return env_miracidia
end



# function to kill humans at age dependent rate


function death_of_human(d)

#= loop through the population, and based on the death rate,
delete individuals from the population  =#

    for i in size(d)[1]:-1:1
        r = rand()
        if ages[i] > 100
            if r < time_step/365
                splice!(d, i)
            end
            else
#= if random number is smaller than the death rate, then delete
    individual from all the arrays  =#
            if r < death_rate[i]
                splice!(d, i)
            end
      end
    end

#=  return the arrays  =#
    return d
end



# function to add a person to the data if a birth occurs

function birth_of_human(d, female_factor, male_factor,
    death_prob_by_age, ages_for_deaths, worm_stages, predis_aggregation,
    mda_adherence, mda_access)


    if rand() > mda_adherence
        adherence = 0
    else
        adherence = 1
    end
    if rand() > mda_access
        access = 0
    else
        access = 1
    end

#  load the gamma distribution for the predispostion distribution
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

    f_worms = fill(0, worm_stages)
    m_worms = fill(0, worm_stages)
    predisp = rand(gamma_pre)[1]
    push!(d, Human(0, rand([0,1]), predisp,
        f_worms, m_worms,
        [], 0, 0, [], [], 1, 0, adherence, access))

    d[end].age_contact_rate = age_contact_rates[1]
    x = d[end].gender
    d[end].predisposition = d[end].predisposition * (1-x) *female_factor + d[end].predisposition * x *female_factor
    d[end].death_age = get_death_age(death_prob_by_age, ages_for_deaths)

    return d
end










# mature the larvae within humans into worms after 35 days


function human_cercariae_maturity(d, time_step)

#=  loop over humans  =#
    for i in 1:length(d)

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_cercariae  =#
        if length(d[i].cercariae) > round(35/time_step, digits = 0)
            females = rand(Binomial(d[i].cercariae[1], 0.5))[1]
            d[i].female_worms[1] = d[i].female_worms[1] + females
            d[i].male_worms[1] = d[i].male_worms[1] + d[i].cercariae[1] - females
            splice!(d[i].cercariae,1)
        end
    end

#= return arrays  =#
    return d
end








# function to administer drug to a specific variable (e.g. female_worms or eggs).
# input the variable, the indices to apply to and the effectiveness of treatment

function administer_drug(d, indices, drug_effectiveness)
    if drug_effectiveness === 1
        @inbounds d[indices] .*= 0
    else
        for i in 1:length(indices)
            @inbounds index = indices[i]
            @inbounds d[index] = rand.(Binomial.(d[index], 1 - drug_effectiveness))
        end
    end
    return d
end







# function for mass drug administration
# currently there is no correlation between individuals chosen each time

function mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender,
    ages, female_worms, male_worms, human_cercariae, eggs,
    treated, mda_round, gender)

#= find index of people with correct ages for the mda treatment =#
    in_gender = in(mda_gender)
    x = findall(((min_age_mda .<= ages .<= max_age_mda) .& in_gender.(gender)))

#=  if this is the first mda round, then treat entirely at random
    find how many people are eligible for the treatment  =#

    #if mda_round == 0
        k = length(x)

#= randomly permute the indices =#
        y = shuffle(x)

#= only take as many as are indicated by the coverage  =#
        y = y[1:trunc(Int, round(k*mda_coverage))]

# update female and male worms, human cercariae and eggs
        female_worms = administer_drug(female_worms, y, mda_effectiveness)
        male_worms = administer_drug(male_worms, y, mda_effectiveness)
        human_cercariae = administer_drug(human_cercariae, y, mda_effectiveness)
        eggs = administer_drug(eggs, y, 1)
 #   else

 #   end

    #println("output = ", female_worms, male_worms, human_cercariae, eggs)
    return female_worms, male_worms, human_cercariae, eggs
    #return x
end




# function to update the mda information


function update_mda(mda_info, mda_round)

    i = min(mda_round + 1, size(mda_info)[1])
    mda_coverage = mda_info[i].coverage
    min_age_mda =  mda_info[i].min_age
    max_age_mda =  mda_info[i].max_age
    mda_effectiveness =  mda_info[i].effectiveness
    mda_gender = mda_info[i].gender
    if mda_round === size(mda_info)[1]
        next_mda_time = Inf
    else
        next_mda_time = mda_info[i].time
    end
    return mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender
end





# function to add vaccination to population

function vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
    vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
    treated, vaccine_duration, vac_status, vaccine_round, gender)

#= find index of people with correct ages for the mda treatment =#
    # x = findall( min_age_vaccine .<= ages .<= max_age_vaccine)

    in_gender = in(vaccine_gender)
    x = findall(((min_age_vaccine .<= ages .<= max_age_vaccine) .& in_gender.(gender)))

#=  if this is the first mda round, then treat entirely at random
    find how many people are eligible for the treatment  =#

    #if mda_round == 0
        k = length(x)

#= randomly permute the indices =#
        y = shuffle(x)

#= only take as many as are indicated by the coverage  =#
        y = y[1:trunc(Int, round(k*vaccine_coverage))]

# update female and male worms, human cercariae and eggs post vaccine
        female_worms = administer_drug(female_worms, y, vaccine_effectiveness)
        male_worms = administer_drug(male_worms, y, vaccine_effectiveness)
        human_cercariae = administer_drug(human_cercariae, y, vaccine_effectiveness)
        eggs = administer_drug(eggs, y, 1)
        vac_status[y] .= vaccine_duration
 #   else

 #   end

    #println("output = ", female_worms, male_worms, human_cercariae, eggs)
    return female_worms, male_worms, human_cercariae, eggs, vac_status
    #return x
end






#= function to update vaccine information =#
function update_vaccine(vaccine_info, vaccine_round)

    i = min(vaccine_round + 1, size(vaccine_info)[1])
    vaccine_coverage = vaccine_info[i].coverage
    min_age_vaccine =  vaccine_info[i].min_age
    max_age_vaccine =  vaccine_info[i].max_age
    vaccine_duration = vaccine_info[i].duration
    vaccine_gender = vaccine_info[i].gender
    if vaccine_round === size(vaccine_info)[1]
        next_vaccine_time = Inf
    else
        next_vaccine_time = vaccine_info[min(vaccine_round + 1, size(vaccine_info)[1])].time
    end
    return vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender

end




function update_env(d, max_fecundity, num_sims, larvae, infective_larvae,
        time_step, contact_rate, r, birth_rate, gamma_pre, female_factor,
        male_factor, age_contact_rates, death_rate_per_time_step)

        sim_time = 0
         update_contact_death_rates = 1
    for j in 1:num_sims
            sim_time += time_step/365
        if sim_time >= update_contact_death_rates
            age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
            update_contact_death_rates += 1
        end
            sim_time += time_step/365
        for i in 1:size(d)[1]
            d[i].age = d[i].age + time_step/365

        end
        # println(1)
        # @time d = larvae_maturity(d, time_step)
        d = human_cercariae_maturity(d, time_step)
        # println(2)
        # @time d = egg_production(d, max_fecundity, r)
        d = egg_production(d, max_fecundity, r)
        # println(3)
        # @time d = worm_maturity(d, worm_stages)
        d = worm_maturity(d, worm_stages)
        # println(4)
        # @time d = vac_decay(d)
        d = vac_decay(d)
        # println(5)
        env_miracidia = miracidia_production(d,  env_miracidia)
        d = death_of_human(d)
        d, infective_larvae, larvae = larvae_uptake(larvae, d, infective_larvae, time_step, contact_rate)



#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_death_rate)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae)

        l = rand(Binomial(size(d)[1], birth_rate))[1]
        if l > 0
            for i in 1:l
                d = birth_of_human(d, gamma_pre, female_factor, male_factor, age_contact_rates,
                death_rate_per_time_step, worm_stages)
            end
        end
    end
    return d
end




function get_prevalences(humans, time)

    pop_burden = [0,0,0]
    sac_burden = [0,0,0]
    adult_burden = [0,0,0]
    pop_prev = 0
    sac_prev = 0
    adult_prev = 0
    sac_pop = 0
    adult_pop = 0
    recorded_eggs =[]
    final_ages = []

    num_humans = length(humans)

    for i in 1:num_humans
        push!(final_ages, humans[i].age);
        final_eggs = humans[i].eggs
        push!(recorded_eggs, final_eggs)
        if final_ages[end] >= 5 && final_ages[end] <= 15
            sac_pop = sac_pop + 1;
        end
        if final_ages[end] > 15
            adult_pop = adult_pop + 1;
        end
        if final_eggs > 16
            pop_burden[3] = pop_burden[3] + 1
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if 5 <= final_ages[end] && final_ages[end] <= 15
                sac_burden[3] = sac_burden[3] + 1
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
            if final_ages[end] > 15
                adult_burden[3] = adult_burden[3] + 1
                adult_burden[2] = adult_burden[2] + 1
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        elseif final_eggs > 4
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if final_ages[end] >= 5 && final_ages[end] <= 15
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
            if final_ages[end] > 15
                adult_burden[2] = adult_burden[2] + 1
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        elseif final_eggs > 0
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if final_ages[end] >= 5 && final_ages[end] <= 15
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev +  1
            end
            if final_ages[end] > 15
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        end
    end

    output = out( round.(100 .*pop_burden./num_humans, digits = 2),
        round.(100 .*sac_burden./sac_pop, digits = 2),
        round.(100 .*adult_burden./adult_pop, digits = 2),
        round.(100 .*pop_prev ./ num_humans, digits = 2),
        round(100 .*sac_prev / sac_pop, digits = 2),
        round(100 .*adult_prev / adult_pop, digits = 2),
        sac_pop, adult_pop,  final_ages, recorded_eggs,
    time)

    return output
end







function update_env_to_equilibrium(num_time_steps,
        humans,
        time_step,
        average_worm_lifespan,
        max_fecundity,
        r,
        worm_stages,
        vaccine_effectiveness,
        density_dependent_fecundity,
        env_miracidia,
        env_cercariae,
        contact_rate,
        env_cercariae_survival_prop,
        env_miracidia_survival_prop,
        record_frequency,
        age_contact_rates,
        human_cercariae_prop)


    sim_time = 0
    record_time = record_frequency
    record = []
    gamma_k = Gamma(0.86,1/0.86)

    d = humans

#=  loop for number of sims  =#

    for j in 1:num_time_steps


        if sim_time >= record_time
            a = get_prevalences(d, sim_time)
            push!(record, a)
            record_time += record_frequency
        end

        sim_time += time_step/365


#=  mature larvae within humans  =#
        d = human_cercariae_maturity(d, time_step)

#=  produce eggs in each human =#
#         eggs = egg_production(eggs, max_fecundity, r, worm_pairs,
#                                total_female_worms, total_male_worms,
#                                density_dependent_fecundity, time_step)

        d = egg_production(d, max_fecundity, r, density_dependent_fecundity, time_step)
#=  mature worms in each human  =#

#         female_worms, male_worms = worm_maturity(female_worms, male_worms,
#                                                  worm_stages, average_worm_lifespan,
#                                                  time_step)

        d = worm_maturity(d, worm_stages, average_worm_lifespan, time_step)


#=  hacth the human eggs into the environment  =#
        # env_miracidia = miracidia_production_by_contact_rate(eggs, env_miracidia, time_step, age_contact_rate)
        env_miracidia = miracidia_production(d, env_miracidia, age_contact_rates)

#=  uptake larvae into humans from the environment  =#
        d, env_cercariae, env_miracidia =
                 cercariae_uptake(d, env_miracidia, env_cercariae, time_step, contact_rate, vaccine_effectiveness)


#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_survival_prop)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_survival_prop)


    end

#=  return the arrays  =#
    return d, env_miracidia, env_cercariae, record
end







function plot_sac_burden_and_sac_high_burden(r, years)
    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac = []
    adult_prev = []

    for i in 1 : length(r)
        push!(times, r[i].time)
        push!(prev, r[i].pop_prev)
        push!(sac_prev, r[i].sac_prev)
        push!(high_burden, r[i].population_burden[3])
        push!(high_burden_sac, r[i].sac_burden[3])
        push!(adult_prev, r[i].adult_prev)
    end

    plot(times, sac_prev, label  = "SAC prevalence", line=(:black, 0.5, 6, :solid))
    plot!(times, high_burden_sac, label  = "SAC high burden", line=(:purple, 0.5, 6))
    plot!(


        size=(800, 600),

        xticks = (0:100:years),
        yticks = 0:10:100,

        ylabel = "Prevalence",
        xlabel = "Year",

        # title  = "The Equation of Time",
        xrotation = rad2deg(pi/3),

        fillrange = 0,
        fillalpha = 0.25,
        fillcolor = :lightgoldenrod,

        background_color = :ivory
        )
    ylims!((0, 100))

end



function get_mean_eggs_age(age_bins, sim_ages, sim_eggs)
    a = []
    for i in 1: length(age_bins)

        if i < length(age_bins)
            x = findall( age_bins[i] .<= sim_ages .<= (age_bins[i+1] - 0.0001))
        else
            x = findall(sim_ages .>= age_bins[i])
        end
        if length(x) > 0
            aa = sim_eggs[x]
            push!(a, mean(aa))
        else
            push!(a, 0)
        end

    end
    return a
end


function get_mean_worms_age(age_bins, sim_ages, male_worms, female_worms)
    a = []
    for i in 1: length(age_bins)

        if i < length(age_bins)
            x = findall( age_bins[i] .<= sim_ages .<= (age_bins[i+1] - 0.0001))
        else
            x = findall(sim_ages .>= age_bins[i])
        end
        if length(x) > 0
            aa = calculate_worm_pairs(female_worms[x], male_worms[x])
            push!(a, mean(aa))
        else
            push!(a, 0)
        end

    end
    return a
end


function repeat_simulations(num_runs, num_sims, N, max_age,
    initial_worms, worm_stages,max_fecundity, infective_larvae, time_step, contact_rate, r,
    birth_rate, gamma_pre, female_factor, male_factor)


    env_miracidia = Int64[]
    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia)
    end
    age_contact_rates = make_age_contact_rate_array(max_age)
    age_death_rate, death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)

    for i in 1:num_runs
        humans = create_population(N, max_age, initial_worms, age_contact_rates,
        death_rate_per_time_step,worm_stages)

        @time humans = update_env(humans, max_fecundity, num_sims, larvae,
        infective_larvae, time_step, contact_rate, r, birth_rate, gamma_pre, female_factor,
        male_factor, age_contact_rates, death_rate_per_time_step)

    j = get_distributions(humans, gamma_k)
    println(j[1])
    end
end



function run_repeated_sims_no_population_change(num_repeats, num_time_steps,
    time_step, average_worm_lifespan,
    max_fecundity, r, worm_stages, predis_aggregation, vaccine_effectiveness,
    density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
    record_frequency, times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, filename,human_cercariae_prop)




    for run in 1:num_repeats


        ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, death_rate_equ, env_miracidia_equ ,
         env_cercariae_equ, adherence_equ, access_equ =
        load_population_from_file(filename, N, true)

        ages , gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,access,
        record =
         update_env_keep_population_same(num_time_steps, copy(ages_equ), copy(human_cercariae_equ), copy(female_worms_equ), copy(male_worms_equ),
                    time_step, average_worm_lifespan,
                    copy(eggs_equ), max_fecundity, r, worm_stages,
                    copy(vac_status_equ), copy(gender_equ), predis_aggregation,
                    copy(predisposition_equ), copy(treated_equ), vaccine_effectiveness,
                    density_dependent_fecundity,
                    copy(vaccinated_equ) , copy(age_contact_rate_equ), copy(death_rate_equ), copy(env_miracidia_equ) ,
                    copy(env_cercariae_equ) , contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                    female_factor, male_factor, contact_rates_by_age,
                    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, copy(adherence_equ), mda_adherence,
                    copy(access_equ), mda_access,
                    record_frequency,human_cercariae_prop);


        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, record, run)

    end
    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev
end








# allow mda and vaccine when updating the model, but whenever there is a death, we instantly add another individual
function update_env_keep_population_same(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
    record_frequency, human_cercariae_prop)

    start_pop = length(ages)
    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = 0
    record = []
    print_time = 0
    if size(mda_info)[1] > 0
        mda_round = 0
        mda_gender = mda_info[1].gender
        mda_coverage = mda_info[1].coverage
        min_age_mda =  mda_info[1].min_age
        max_age_mda =  mda_info[1].max_age
        mda_effectiveness =  mda_info[1].effectiveness
        next_mda_time = mda_info[1].time
    else
        next_mda_time = Inf
    end


    if size(vaccine_info)[1] > 0
        vaccine_round = 0
        vaccine_coverage = vaccine_info[1].coverage
        vaccine_gender = vaccine_info[1].gender
        min_age_vaccine =  vaccine_info[1].min_age
        max_age_vaccine =  vaccine_info[1].max_age
        next_vaccine_time = vaccine_info[1].time
        vaccine_duration = vaccine_info[1].duration
    else
        next_vaccine_time = Inf
    end
#=  loop for number of sims  =#

    for j in 1:num_time_steps

#= update contact and death rates every year =#
        if sim_time >= update_contact_death_rates
            death_rate = update_death_rate(ages, death_rate, death_rate_per_time_step)
            age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences(d, sim_time)
            push!(record, a)
            record_time += record_frequency
        end

        sim_time += time_step/365


#=  mature larvae within humans  =#
   human_cercariae, female_worms, male_worms =
            human_cercariae_maturity(human_cercariae, female_worms, male_worms, time_step)


#=  calculate the number of worm pairs in each human  =#
        worm_pairs = calculate_worm_pairs(female_worms, male_worms)

#=  calculate the total number of worms in each human  =#
       total_female_worms, total_male_worms =
        calculate_total_worms(female_worms, male_worms)

#=  produce eggs in each human =#
        eggs = egg_production(eggs, max_fecundity, r, worm_pairs,
                               total_female_worms, total_male_worms,
                               density_dependent_fecundity, time_step)

#=  mature worms in each human  =#
        female_worms, male_worms = worm_maturity(female_worms, male_worms,
                                                 worm_stages, average_worm_lifespan,
                                                 time_step)

 #=  reduce the vaccination status by the time step  =#
        vac_status = vac_status .- time_step/365

#=  hacth the human eggs into the environment  =#
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step, age_contact_rate)

    #=  update population due to death  =#
       ages , gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate, adherence, access =
             death_of_human(ages, gender, predisposition,  human_cercariae, eggs,
                                         vac_status, treated, female_worms, male_worms,
                                         vaccinated, age_contact_rate, death_rate, time_step,
                                         adherence, access)

        if length(ages) < start_pop
            num_births = start_pop - length(ages)
            for i in 1:num_births
                ages , gender, predisposition,  human_cercariae, eggs, vac_status,
                 treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate, adherence, access =
                     birth_of_human(ages , gender, predisposition,  human_cercariae, eggs, vac_status,
                                        treated, female_worms, male_worms,vaccinated, age_contact_rate,
                                        death_rate, female_factor, male_factor, contact_rates_by_age,
                                        death_rate_per_time_step, worm_stages, predis_aggregation, adherence,
                                        mda_adherence, access, mda_access)
            end
        end
#=  uptake larvae into humans from the environment  =#
        env_cercariae, human_cercariae, env_miracidia =
                 cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            female_worms, male_worms, human_cercariae, eggs =
            mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender,
                     ages, female_worms, male_worms, human_cercariae, eggs,
                     treated, mda_round, gender, adherence, access)

#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            female_worms, male_worms, human_cercariae, eggs, vac_status =
                    vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                    vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
                    treated, vaccine_duration, vac_status, vaccine_round, gender, adherence, access)

#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                update_vaccine(vaccine_info, vaccine_round)
        end

#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_survival_prop)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_survival_prop, time_step)


    end

#=  return the arrays  =#
    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,access,
    record
end




using Distributions
using Random
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
contact_rate = 0.007
max_fecundity = 0.34  # From "The design of schistosomiasis monitoring and evaluation programmes:
#The importance of collecting adu lt data to inform treatment strategies for Schistosoma mansoni"
density_dependent_fecundity = 0.0007
r = 0.03 # aggregation parameter for negative binomial for egg production
num_time_steps = trunc(Int, 365*50/ time_step)

#const num_time_steps =
birth_rate = 20*time_step/(1000*365)

average_worm_lifespan = 5.7 # years
predis_aggregation = 0.24
env_cercariae_death_rate = 0.09 * time_step #= life span of cercariae in the environment is short 8-20 hrs
according to "Studies of the Transmission Dynamics, Mathematical Model Development and the Control of Schistosome
Parasites by Mass Drug Administration in Human Communities"  =#
env_miracidia_death_rate = 0 * time_step
mda_coverage = 0.8 # proportion of target age group reached by mda
mda_round = 0
num_runs = 1
env_cercariae = 0
predis_aggregation = 0.24
gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

vaccine_effectiveness = 0.86

env_cercariae_survival_prop = 1
env_miracidia_survival_prop = 1

record_frequency = 1/24
human_cercariae_prop = 1
#=
    create array with the number of initial larvae which have been living for as long
    as initial_larvae_days
=#
env_miracidia = Int64[]
for i in 1 : initial_miracidia_days
    push!(env_miracidia, initial_miracidia)
end

death_prob_by_age = [0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                            0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1]

    ages_for_deaths = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                          65, 70, 75, 80, 85, 90, 95, 100, 110];



#=
    Create population. Initialize with random ages,
    random gender (0 = female, 1 = male), predisposition
=#


mda_adherence = 0.8
mda_access = 0.8
scenario = "high adult"
spec_ages = 7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
            2636, 1970, 1468, 1166, 943, 718, 455, 244
ages_per_index = 5
ages = specified_age_distribution(N, spec_ages, ages_per_index)
age_contact_rates = make_age_contact_rate_array(max_age, scenario, [], []);

N = 2000
humans, env_miracidia = create_population_specified_ages(N, initial_worms, age_contact_rates,
        worm_stages, female_factor, male_factor,initial_miracidia,
        initial_miracidia_days, predis_aggregation, time_step,
        spec_ages, ages_per_index,
        mda_adherence, mda_access);

        N = 1000



        contact_rate = 0.035
        r = 0.03 # aggregation parameter for negative binomial for egg production
        human_cercariae_prop = 1
        env_cercariae_survival_prop = 1/10
        env_miracidia_survival_prop = 1/12

        num_years_equ = 500
        num_time_steps = trunc(Int, 365*num_years_equ/ time_step)

        humans, env_miracidia = create_population_specified_ages(N, initial_worms, age_contact_rates,
                worm_stages, female_factor, male_factor,initial_miracidia,
                initial_miracidia_days, predis_aggregation, time_step,
                spec_ages, ages_per_index,
                mda_adherence, mda_access);

        @time humans, env_miracidia, env_cercariae, record = update_env_to_equilibrium(num_time_steps,
                humans,
                time_step,
                average_worm_lifespan,
                max_fecundity,
                r,
                worm_stages,
                vaccine_effectiveness,
                density_dependent_fecundity,
                env_miracidia,
                env_cercariae,
                contact_rate,
                env_cercariae_survival_prop,
                env_miracidia_survival_prop,
                record_frequency,
                age_contact_rates,
                human_cercariae_prop);


        plot_sac_burden_and_sac_high_burden(record, num_years_equ)
