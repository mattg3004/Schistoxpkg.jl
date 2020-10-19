using Distributions
using Random
using JLD


mutable struct out
    population_burden
    sac_burden
    adult_burden
    pop_prev
    sac_prev
    adult_prev
    sac_pop
    adult_pop
    final_ages
    recorded_eggs
    time
end


mutable struct mda_information
    coverage
    min_age
    max_age
    gender
    effectiveness
    time
end


mutable struct vaccine_information
    coverage
    min_age
    max_age
    gender
    duration
    time
end



mutable struct Human
    age::Float64
    death_age::Float64
    gender::Int64
    predisposition::Float64
    female_worms::Array{Int64}
    male_worms::Array{Int64}
    eggs::Int64
    vac_status::Int64
    age_contact_rate::Float64
    adherence::Int64
    access::Int64
    community::Int64
    relative_contact_rate::Float64
    uptake_rate::Float64
    acquired_immunity::Float64 # level of acquired immunity
    total_worms::Int64 # total number of worms over lifetime
    larvae::Array{Int64}
    Human() = Human(0.0,100, 0,1.0,[],[],0,0,0.0,0,0,0,0, 0, 0, 0,[])
    Human(age,death_age, gender, predisposition, female_worms, male_worms, eggs, vac_status,
        age_contact_rate, adherence, access, community, relative_contact_rate, uptake_rate,
        acquired_immunity, total_worms, larvae) =
    new(age, death_age, gender, predisposition, female_worms, male_worms, eggs, vac_status,
        age_contact_rate, adherence, access, community, relative_contact_rate, uptake_rate,
        acquired_immunity, total_worms, larvae)
end



mutable struct Parameters
    N::Int64
    time_step::Float64
    N_communities::Int64
    community_probs::Array{Float64}
    community_contact_rate::Array{Float64}
    density_dependent_fecundity::Float64
    average_worm_lifespan::Float64
    max_age::Float64
    initial_worms::Int64
    initial_miracidia::Int64
    initial_miracidia_days::Int64
    init_env_cercariae::Int64
    worm_stages::Int64
    contact_rate::Float64
    max_fecundity::Float64
    age_contact_rates::Array{Float64}
    ages_for_contacts::Array{Int64}
    contact_rate_by_age_array::Array{Float64}
    mda_adherence::Float64
    mda_access::Float64
    female_factor::Float64
    male_factor::Float64
    miracidia_maturity::Int64
    birth_rate::Float64
    human_cercariae_prop::Float64
    predis_aggregation::Float64
    cercariae_survival::Float64
    miracidia_survival::Float64
    death_prob_by_age::Array{Float64}
    ages_for_death::Array{Float64}
    r::Float64
    vaccine_effectiveness::Float64
    drug_effectiveness::Float64
    spec_ages::Array{Float64}
    ages_per_index::Int64
    record_frequency::Float64
    use_kato_katz::Int64 # if 0, then don't use KK, if 1, use KK
    kato_katz_par::Float64
    heavy_burden_threshold::Int64
    rate_acquired_immunity::Float64
    M0::Float64 # mean worm burden
    human_larvae_maturity_time::Int64

    # constructors with default values
    Parameters() =
          Parameters(1000,1, 1, [1.0], [1.0],
          0.0007, 5.7, 100, 0,
          1000, 2, 1000,
          1,  0.03, 0.34,[0.032,0.61, 1,0.06],[4,9,15,100],fill(0,101),
          # gamma_r::Distributions.Gamma{Float64}     # distribution of lifelong predisposition factor
          # gamma_k::Distributions.Gamma{Float64}     # distribution for kato_katz test
          1.0, 1.0, 1.0,
          1.0, 24, 1.0, 1.0, 0.24, 0.5, 0.5,
        [0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                     0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1],
        [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                   65, 70, 75, 80, 85, 90, 95, 100, 110], 0.03, 1,1,
    [7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
              2636, 1970, 1468, 1166, 943, 718, 455, 244],5,1/24, 0, 0.87, 16, 0, 15, 0)

    Parameters(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days, init_env_cercariae,
        worm_stages, contact_rate, max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival, miracidia_survival,
        death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index, record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
        rate_acquired_immunity, M0, human_larvae_maturity_time) =
    new(N, time_step, N_communities, community_probs, community_contact_rate,
        density_dependent_fecundity, average_worm_lifespan,
        max_age, initial_worms, initial_miracidia, initial_miracidia_days,
        init_env_cercariae, worm_stages, contact_rate, max_fecundity, age_contact_rates,
        ages_for_contacts, contact_rate_by_age_array, mda_adherence, mda_access,  female_factor, male_factor, miracidia_maturity,
        birth_rate, human_cercariae_prop, predis_aggregation,cercariae_survival, miracidia_survival,
                death_prob_by_age, ages_for_death, r, vaccine_effectiveness, drug_effectiveness,
        spec_ages, ages_per_index,record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
        rate_acquired_immunity, M0, human_larvae_maturity_time)
end



mutable struct Environment
     miracidia::Array{Int64}
     cercariae::Int64
     humans::Array{Human}          # container for population
     time::Int64                     # in years
     record::Array{out}
     # Constructors
     Environment() = Environment([],0,[],0, [])
     Environment(miracidia, cercariae, humans, time, record) = new(miracidia, cercariae, humans, time, record)
end




# create the age specific contact settings given the scenario
function create_contact_settings(scenario)
    if scenario == "low adult"
        age_contact_rates = [0.01, 1.2, 1, 0.02]
    elseif scenario == "moderate adult"
        age_contact_rates = [0.032, 0.61, 1, 0.06]
    elseif scenario == "high adult"
        age_contact_rates = [0.01, 0.61, 1, 0.12]
    end
    return age_contact_rates
end



# function to get age dependent contact rate.
# the contact rates are taken from the
# "What is required in terms of mass drug administration to interrupt the transmission
#     of schistosome parasites in regions of endemic infection?" paper
# at some point we may change this to be an input from a file instead

function make_age_contact_rate_array(pars, scenario, input_ages, input_contact_rates)
    if pars.max_age < 60
        error("max_age must be greater than 60")
    else

    if length(input_ages) == 0
        contact_settings = create_contact_settings(scenario)
        x =[fill(contact_settings[4], trunc(Int,(pars.max_age+1)))]
    # initialize an array with the same value for contact rate across all ages
        pars.contact_rate_by_age_array = fill(contact_settings[4], trunc(Int,(pars.max_age+1)))
#         pars.contact_rate_by_age_array = pars.contact_rate_by_age_array[1]

    # then edit the entries for different ages according to values
        for i in 1:5
            pars.contact_rate_by_age_array[i] = contact_settings[1]
        end
        if scenario == "high adult"
            for i in 6:12
                pars.contact_rate_by_age_array[i] = contact_settings[2]
            end
            for i in 13:21
                pars.contact_rate_by_age_array[i] = contact_settings[3]
            end
        else
            for i in 6:10
                pars.contact_rate_by_age_array[i] = contact_settings[2]
            end
            for i in 11:16
                pars.contact_rate_by_age_array[i] = contact_settings[3]
            end
        end

        else
            for i in 1:length(pars.contact_rate_by_age_array)
                pars.contact_rate_by_age_array[i] = input_contact_rates[end]
            end

            for i in 1 : length(input_contact_rates)
                if i == 1
                    for j in 1:trunc(Int,(input_ages[i] + 1))
                        pars.contact_rate_by_age_array[j] = input_contact_rates[i]
                    end
                else
                    for j in trunc(Int,(input_ages[i-1]+2)):trunc(Int,(input_ages[i] + 1))
                        pars.contact_rate_by_age_array[j] = input_contact_rates[i]
                    end
                end
            end

        end

        return pars
    end


# return contact rates for age
end



# function to generate an age for death of an individual
function get_death_age(pars)
        age = 0
        k = 1
        p = pars.death_prob_by_age[k]
        goal_age = pars.ages_for_death[k]
        x = rand()
        while x > p
            age += 1
            if age >= goal_age
                k += 1
                p = pars.death_prob_by_age[k]
                goal_age = pars.ages_for_death[k]
            end
            x = rand()
        end
        if k == 1
            min_age = 0
        else
            min_age = pars.ages_for_death[k-1]
        end

        death_age = min_age + rand() * (goal_age - min_age)

    return death_age
end



# function to age population and generating death ages
function generate_ages_and_deaths(num_steps, humans, pars)
     for i in 1:num_steps
        x = []
        for i in 1:length(humans)
            if( humans[i].age - humans[i].death_age > 0)
                push!(x, i)
            end
            humans[i].age += pars.time_step/365
        end
        k = length(x)

        if k > 0
            for j in k:-1:1
                splice!(humans, x[j])
                humans = birth_of_human(humans, pars)
            end
        end
    end
    return humans
end




#=
    define a function which will create the initial population.
    This randomly chooses age, and gender.
    predisposition is taken to be gamma distributed.
    There is also a male and female adjustment to predisposition
    adjusting for gender specific behaviour
=#
"""
    create_population

This will create the initial human population with randomly chosen age, and gender.
Predisposition is taken to be gamma distributed
There is also a male and female adjustment to predisposition adjusting for gender specific behaviour
In addition to this, it will create the initial miracidia environment vector
"""
function create_population(pars)


     if length(pars.community_probs) != pars.N_communities
        error("must provide probabilities for membership of each community")
    else
        community_selection = 1
    end
    if pars.N_communities > 1
        community_selection = cumsum(pars.community_probs)/sum(pars.community_probs)
    end

    env = Environment()
    humans = env.humans

#=  initialize the Gamma distribution for predisposition selection  =#
    pre = Gamma(pars.predis_aggregation, 1/pars.predis_aggregation)

#= select all predispositions =#
    predisposition = rand(pre, pars.N)

#=  initialize and fill the environmental variable  =#

    miracidia = env.miracidia
    for i in 1 : pars.initial_miracidia_days
        push!(miracidia, pars.initial_miracidia)
    end
    cercariae = env.cercariae
    cercariae = pars.init_env_cercariae;
    for i in 1:pars.N

        f_worms = fill(0, pars.worm_stages)
        f_worms[1] = trunc(Int,round(rand()*pars.initial_worms))
        m_worms = fill(0, pars.worm_stages)
        m_worms[1] = trunc(Int,round(rand()*pars.initial_worms))

        if rand() > pars.mda_adherence
            adherence = 0
        else
            adherence = 1
        end

        if rand() > pars.mda_access
            access = 0
        else
            access = 1
        end

        community = findall(community_selection .> rand())[1]

# Human(age,gender, predisposition, female_worms, male_worms, eggs, vac_status, age_contact_rate, adherence, access, community)
        death_age = get_death_age(pars)
        push!(humans, Human(pars.max_age*rand(), death_age, rand([0,1]), predisposition[i],
        f_worms, m_worms,
        0, 0, 0, adherence, access, community, 0, 0,0,0, []))

        age = trunc(Int, humans[end].age)

        contact_rate_age = findall(pars.ages_for_contacts .> age)[1]


        humans[end].age_contact_rate = pars.age_contact_rates[contact_rate_age]
        humans[end].relative_contact_rate = humans[end].age_contact_rate *  pars.community_contact_rate[humans[end].community]/
        (maximum(pars.community_contact_rate) * maximum(pars.contact_rate_by_age_array))


        if humans[end].gender == 0
            humans[end].predisposition = humans[end].predisposition * pars.female_factor
        else
            humans[end].predisposition = humans[end].predisposition * pars.male_factor
        end

        humans[end].uptake_rate = humans[end].predisposition * pars.contact_rate * humans[end].age_contact_rate *
                                    pars.time_step * pars.community_contact_rate[community]
    end
    return humans, miracidia, cercariae

end





#=
    define a function which will create the initial population.
    This randomly chooses age, and gender.
    predisposition is taken to be gamma distributed.
    There is also a male and female adjustment to predisposition
    adjusting for gender specific behaviour
=#
"""
    create_population_specified_ages(pars)

This will create the initial human population with an age distribution
specified by the spec_ages variable
Predisposition is taken to be gamma distributed. There is also a male and female
adjustment to predisposition adjusting for gender specific behaviour
In addition to this, it will create the initial miracidia environment vector
"""
function create_population_specified_ages(pars)


     if length(pars.community_probs) != pars.N_communities
        error("must provide probabilities for membership of each community")
    else
        community_selection = 1
    end
    if pars.N_communities > 1
        community_selection = cumsum(pars.community_probs)/sum(pars.community_probs)
    end
    ages = specified_age_distribution(pars)
    env = Environment()
    humans = env.humans

#=  initialize the Gamma distribution for predisposition selection  =#
    pre = Gamma(pars.predis_aggregation, 1/pars.predis_aggregation)

#= select all predispositions =#
    predisposition = rand(pre, pars.N)

#=  initialize and fill the environmental variable  =#

    miracidia = env.miracidia
    for i in 1 : pars.initial_miracidia_days
        push!(miracidia, pars.initial_miracidia)
    end
    cercariae = env.cercariae
    cercariae = pars.init_env_cercariae;
    for i in 1:pars.N

        f_worms = fill(0, pars.worm_stages)
        f_worms[1] = trunc(Int,round(rand()*pars.initial_worms))
        m_worms = fill(0, pars.worm_stages)
        m_worms[1] = trunc(Int,round(rand()*pars.initial_worms))

        if rand() > pars.mda_adherence
            adherence = 0
        else
            adherence = 1
        end

        if rand() > pars.mda_access
            access = 0
        else
            access = 1
        end

        community = findall(community_selection .> rand())[1]

# Human(age,gender, predisposition, female_worms, male_worms, eggs, vac_status, age_contact_rate, adherence, access, community)
        death_age = get_death_age(pars)
        push!(humans, Human(ages[i]+rand(), death_age, rand([0,1]), predisposition[i],
        f_worms, m_worms,
        0, 0, 0, adherence, access, community, 0, 0,0,0, []))

        age = trunc(Int, humans[end].age)

        contact_rate_age = findall(pars.ages_for_contacts .> age)[1]


        humans[end].age_contact_rate = pars.age_contact_rates[contact_rate_age]
        humans[end].relative_contact_rate = humans[end].age_contact_rate *  pars.community_contact_rate[humans[end].community]/
        (maximum(pars.community_contact_rate) * maximum(pars.contact_rate_by_age_array))


        if humans[end].gender == 0
            humans[end].predisposition = humans[end].predisposition * pars.female_factor
        else
            humans[end].predisposition = humans[end].predisposition * pars.male_factor
        end

        humans[end].uptake_rate = humans[end].predisposition * pars.contact_rate * humans[end].age_contact_rate *
                                    pars.time_step * pars.community_contact_rate[community]
    end
    return humans, miracidia, cercariae

end







"""
    update_contact_rate(humans, pars)

function to update the contact rate of individuals in the population. This is necessary
    as over time when people age, they will move through different age groups which have
    different contact rates
"""
function update_contact_rate(humans,  pars)
    @inbounds for h in humans
        age = min(length(pars.contact_rate_by_age_array) - 1, (trunc(Int, h.age)))
        h.age_contact_rate =  pars.contact_rate_by_age_array[age+1]
        h.relative_contact_rate = h.age_contact_rate *  pars.community_contact_rate[h.community]/
        (maximum(pars.community_contact_rate) * maximum(pars.contact_rate_by_age_array))

        h.uptake_rate = h.predisposition * pars.contact_rate * h.age_contact_rate *
                                    pars.time_step * pars.community_contact_rate[h.community]
    end
    return humans
end




#=
    humans uptake larvae based on their predisposition, age_dependent contact rate
    and the number of larvae in the environment
    input d here is the whole humans Array
=#


"""
    cercariae_uptake(humans, cercariae, miracidia, pars)

uptake cercariae into humans, whilst updating cercariae with miracidia.
Uptaken cercariae immediately become worms in this formulation
"""
function cercariae_uptake!(humans, cercariae, miracidia, pars)

    humans = shuffle(humans)
    k = length(humans)
#= assign larvae which have been in the environment for 40 days to become infective.
    then delete those larvae from the environmental larvae =#

    if  length(miracidia) > (pars.miracidia_maturity / pars.time_step)
        cercariae += (miracidia[1] * pars.human_cercariae_prop)
        splice!(miracidia, 1)
    end

    @inbounds for h in humans
#= if there are still infective larvae in the environment,
we will uptake from choose from Poisson
distribution. otherwise, just uptake 0. =#
     #   if cercariae > 0
# calculate the rate of the poisson distribution
        pois_rate  = max(h.uptake_rate *  cercariae / k, 0)

        # reduce the rate according to the effectiveness of the vaccine (if any is given)
        pois_rate = pois_rate * (1 - (h.vac_status > 0) * pars.vaccine_effectiveness) * (1 - h.acquired_immunity)

        # choose from the Poisson distribution
        uptake = rand(Poisson(pois_rate))
            #println(uptake)
        n = rand(Binomial(uptake, 0.5))

        h.female_worms[1] += n
        h.male_worms[1] += uptake - n
        h.total_worms += uptake

        # # reduce the infective larvae by the number of larvae uptaken
        cercariae -= uptake
        cercariae = max(0, cercariae)
    end
    return humans, cercariae, miracidia
end



"""
    cercariae_uptake_human_larvae!(humans, cercariae, miracidia, pars)

uptake cercariae into humans, whilst updating cercariae with miracidia.
Uptaken cercariae become larvae within humans, rather than immmediately into worms with this function.
"""
function cercariae_uptake_with_human_larvae!(humans, cercariae, miracidia, pars)

    humans = shuffle(humans)
    k = length(humans)
#= assign larvae which have been in the environment for 40 days to become infective.
    then delete those larvae from the environmental larvae =#

    if  length(miracidia) > (pars.miracidia_maturity / pars.time_step)
        cercariae += (miracidia[1] * pars.human_cercariae_prop)
        splice!(miracidia, 1)
    end

    @inbounds for h in humans
#= if there are still infective larvae in the environment,
we will uptake from choose from Poisson
distribution. otherwise, just uptake 0. =#
     #   if cercariae > 0
# calculate the rate of the poisson distribution
        pois_rate  = max(h.uptake_rate * (1-pars.rate_acquired_immunity * h.total_worms) *  cercariae / k, 0)

        # reduce the rate according to the effectiveness of the vaccine (if any is given)
        pois_rate = pois_rate * (1 - (h.vac_status > 0) * pars.vaccine_effectiveness) * (1 - h.acquired_immunity)

        # choose from the Poisson distribution
        uptake = rand(Poisson(pois_rate))
            #println(uptake)
        push!(h.larvae, uptake)

        # # reduce the infective larvae by the number of larvae uptaken
        cercariae -= uptake
        cercariae = max(0, cercariae)
    end
    return humans, cercariae, miracidia
end



"""
    human_larvae_maturity(humans, pars)

This will mature the human larvae into worms after a chosen number of days, which is specified
by the human_larvae_maturity_time parameter in the pars struct
"""
function human_larvae_maturity(humans, pars)

#=  loop over humans  =#
    @inbounds for h in humans

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_cercariae  =#
        if length(h.larvae) > round(pars.human_larvae_maturity_time/time_step, digits = 0)
            females = rand(Binomial(h.larvae[1], 0.5))
            h.female_worms[1] += females
            h.male_worms[1] += h.larvae[1] - females
            splice!(h.larvae, 1)
        end
    end

#= return arrays  =#
    return humans
end




# function to kill miracidia in the environment
"""
    miracidia_death!(miracidia, pars)

Kill a chosen proportion of miracidia in the environment governed by the
miracidia_survival parameter in the pars struct
"""
function miracidia_death!(miracidia, pars)
    #= as emiracidia is an array, we need to use . syntax to apply
    the functions to each element in the array =#
    # return rand.(Binomial.(env_miracidia, 1 - env_miracidia_survival_prop))
    if pars.miracidia_survival <= 0
        error("miracidia_survival_prop must be bigger than 0")
    else
        miracidia[end] = trunc(Int, round(miracidia[end] * pars.miracidia_survival, digits = 0))
    end
    return miracidia
end




# function to kill cercariae in the environment
"""
    cercariae_death!(miracidia, pars)

Kill a chosen proportion of cercariae in the environment governed by the
cercariae_survival parameter in the pars struct
"""
function cercariae_death!(cercariae, pars)
    # updated_cercariae = 0
    # for i in 1:time_step
    #     updated_cercariae += (env_cercariae/time_step) * (1/(2^i))
    # end
    if pars.cercariae_survival <= 0
        error("cercariae_survival_prop must be bigger than 0")
    else
        updated_cercariae = trunc(Int, round(cercariae * pars.cercariae_survival, digits= 0))
    end
    return updated_cercariae
end






# Worms die have a specified mean life span, and hence a rate of deaths per day .
# the number of deaths, or maturing from each stage is dependent
# on the number of worm stages and rate of aging through stages and dying.
# This p is multiplied by the time scale of the simulation, so if 2 days
# pass between consecutive time points, twice as many worms age and die

# for human in humans
#     println(human.eggs)
# end
"""
    worm_maturity!(humans, pars)

function to kill worms within human hosts, and if there is more than one stage for worm life,
to update how many worms are in each stage
"""
function worm_maturity!(humans, pars)

    # probability of aging out of category/ dying
    p = pars.time_step * pars.worm_stages/ (365 * pars.average_worm_lifespan)

    @inbounds for h in humans

        # kill appropriate number of worms in the final stage
         n = h.female_worms[pars.worm_stages]
        h.female_worms[pars.worm_stages] = rand(Binomial(n, 1-p))

        n = h.male_worms[pars.worm_stages]
        h.male_worms[pars.worm_stages] = rand(Binomial(n, 1-p))

#=
        for aging worms, we do this in reverse order, which ensures the
        correct order of aging is respected
=#

        @inbounds for j in (pars.worm_stages-1):-1:1
#=   choose the number of male and female worms to age from one stage to the next   =#
            aging_females = rand(Binomial(h.female_worms[j], p))
            aging_males = rand(Binomial(h.male_worms[j], p))

#=   add and subtract the number of worms from the appropriate categories   =#
            h.female_worms[j+1] += aging_females
            h.female_worms[j] -= aging_females
            h.male_worms[j+1] += aging_males
            h.male_worms[j] -= aging_males
        end
    end
    return humans

end


# function calculate_worm_pairs(humans)
#     return min(sum(humans.female_worms), sum(humans.male_worms))
# end

"""
    calculate_worm_pairs(female_worms, male_worms)

calculate how many pairs of worms there are in each human host
"""
function calculate_worm_pairs(female_worms, male_worms)
    worm_pairs = Int64[]
    for i in 1:length(female_worms)
        push!(worm_pairs, min(sum(female_worms[i]), sum(male_worms[i])))
    end
    return(worm_pairs)
    #return min.(sum.(female_worms), sum.(male_worms))
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

"""
    egg_production!(humans, pars)

function to produce eggs for individuals, dependent on how many worms they have
        and the max fecundity and density dependent fecundity of the population
"""
function egg_production!(humans, pars)
    male_worms = (p->p.male_worms).(humans)
    female_worms = (p->p.female_worms).(humans)
    worm_pairs = calculate_worm_pairs(female_worms, male_worms)

    for i in 1 : length(humans)
#        wp = calculate_worm_pairs(humans[i])
        wp = worm_pairs[i]
        wp = max(wp,0.000000001)
#= if we have a positive number of worms, then make calculation,
        otherwise the number of eggs is trivially 0 =#
#         if worm_pairs > 0

                     #println(female_worms[i][1])
# calculate the mean number of eggs we would expect
    #    mean_eggs = pars.max_fecundity * wp *
    #            exp(- pars.density_dependent_fecundity  * female_worms[i][1])
    mean_eggs = pars.max_fecundity * wp *
                exp(- pars.density_dependent_fecundity  * wp)


        # mean_eggs = 0.5*(pars.max_fecundity * M0)/(M0 + M) * (1-exp(-M0*log(2)*M)) * M
# calculate the number of successes
    #    NB_r = pars.r * wp

# calculate the probability of a success
    #    p = NB_r/(NB_r+mean_eggs)

# choose from NB
         #eggs = rand(NegativeBinomial(NB_r,p))
        eggs = rand(Poisson(mean_eggs))
            # println(eggs)
        humans[i].eggs = eggs
    end
    return humans
end



"""
    egg_production_increasing!(humans, pars)

function to produce eggs for individuals, dependent on how many worms they have
        and the max fecundity and density dependent fecundity of the population
"""
function egg_production_increasing!(humans, pars)
    male_worms = (p->p.male_worms).(humans)
    female_worms = (p->p.female_worms).(humans)
    worms = sum.(male_worms) + sum.(female_worms)
    for i in 1 : length(humans)
#        wp = calculate_worm_pairs(humans[i])
        M = worms[i]
#= if we have a positive number of worms, then make calculation,
        otherwise the number of eggs is trivially 0 =#
#         if worm_pairs > 0

                     #println(female_worms[i][1])
# calculate the mean number of eggs we would expect


         mean_eggs = 0.5*(pars.max_fecundity * pars.M0)/(pars.M0 + M) * (1-exp(-pars.M0*log(2)*M)) * M
# calculate the number of successes


# calculate the probability of a success


# choose from NB
         #eggs = rand(NegativeBinomial(NB_r,p))
        eggs = rand(Poisson(mean_eggs))
            # println(eggs)
        humans[i].eggs = eggs
    end
    return humans
end





# function
"""
    miracidia_production!(humans)

release eggs from individuals into the environment as miracidia. Release is relative to the contact rate with
the environment for each individual.
"""
function miracidia_production!(humans)
    released_eggs = 0
    @inbounds for h in humans
        released_eggs += h.eggs * h.relative_contact_rate
    end
    return (round(released_eggs))
end



"""
    death_of_human(humans)

remove individuals from the population whose age is greater than their death age
"""
function death_of_human(humans)

    for i in length(humans):-1:1
        if humans[i].age > humans[i].death_age
            splice!(humans, i)
        end
    end
    return humans
end



"""
    birth_of_human(humans, pars)

add an individual to the population
"""
function birth_of_human(humans, pars)
    if length(pars.community_probs) != pars.N_communities
        error("must provide probabilities for membership of each community")
    else
        community_selection = 1
    end
    if pars.N_communities > 1
        community_selection = cumsum(pars.community_probs)/sum(pars.community_probs)
    end

    pre = Gamma(pars.predis_aggregation, 1/pars.predis_aggregation)
    predisp = rand(pre)[1]
    f_worms = fill(0, pars.worm_stages)
    m_worms = fill(0, pars.worm_stages)


    if rand() > pars.mda_adherence
        adherence = 0
    else
        adherence = 1
    end

    if rand() > pars.mda_access
        access = 0
    else
        access = 1
    end

    community = findall(community_selection .> rand())[1]

# Human(age,gender, predisposition, female_worms, male_worms, eggs, vac_status, age_contact_rate, adherence, access, community)
    death_age = get_death_age(pars)


    push!(humans, Human(0, death_age, rand([0,1]), predisp,
        f_worms, m_worms,
        0, 0, 0, adherence, access, community,0,0,0,0, []))


    humans[end].age_contact_rate = pars.age_contact_rates[1]



    humans[end].relative_contact_rate = humans[end].age_contact_rate *  pars.community_contact_rate[humans[end].community]/
        (maximum(pars.community_contact_rate) * maximum(pars.contact_rate_by_age_array))


    if humans[end].gender == 0
        humans[end].predisposition = humans[end].predisposition * pars.female_factor
    else
        humans[end].predisposition = humans[end].predisposition * pars.male_factor
    end

    humans[end].uptake_rate = humans[end].predisposition * pars.contact_rate * humans[end].age_contact_rate *
                                    pars.time_step * pars.community_contact_rate[community]

    return humans
end



# function to administer drug to a specific variable (e.g. female_worms or eggs).
# input the variable, the indices to apply to and the effectiveness of treatment


"""
    administer_drug(humans, indices, drug_effectiveness)

administer mda drugs to chosen individuals in the population. If they adhere to the drugs, then they
reduce male and female worms with a given efficacy alongside removing eggs
"""
function administer_drug(humans, indices, drug_effectiveness)

    @inbounds for i in 1:length(indices)
         index = indices[i]
        p = 1 - (drug_effectiveness * humans[index].adherence)
        humans[index].female_worms = rand.(Binomial.(humans[index].female_worms,p))
        humans[index].male_worms = rand.(Binomial.(humans[index].male_worms, p))
#         @inbounds humans[index].human_cercariae = rand.(Binomial.(humans[index].human_cercariae,p))
        humans[index].eggs = 0
    end
    return humans
end


# function to administer drug to a specific variable (e.g. female_worms or eggs).
# input the variable, the indices to apply to and the effectiveness of treatment
"""
    administer_vaccine(humans, indices, vaccine_effectiveness, vaccine_duration)

administer vaccine to chosen individuals in the population.
reduce male and female worms with a given efficacy alongside removing eggs and
adding to their vaccine status signifying that they will have increased immunity for a chosen period of time
"""
function administer_vaccine(humans, indices, vaccine_effectiveness, vaccine_duration)

    @inbounds for i in 1:length(indices)
         index = indices[i]
        p = 1 - (vaccine_effectiveness)
        humans[index].female_worms = rand.(Binomial.(humans[index].female_worms,p))
        humans[index].male_worms = rand.(Binomial.(humans[index].male_worms, p))
        humans[index].eggs = 0
        humans[index].vac_status = vaccine_duration * humans[index].adherence
    end
    return humans
end


# function for mass drug administration
# currently there is no correlation between individuals chosen each time
"""
    mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)

administer mda in the population. This includes choosing individuals between specified ages,
having a certain level of coverage and taking access and adherence into consideration
"""
function mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)

    ages = (p->p.age).(humans)
    gender = (p->p.gender).(humans)
    access = (p->p.access).(humans)

    #= find index of people with correct ages for the mda treatment =#
    in_gender = in(mda_gender)
    x = findall(((min_age_mda .<= ages .<= max_age_mda) .& in_gender.(gender) .& (access .== 1)))


    k = length(x)


#= randomly permute the indices =#
    y = shuffle(x)

#= only take as many as are indicated by the coverage  =#
    y = y[1:trunc(Int, round(k*mda_coverage))]

# update female and male worms, human cercariae and eggs

    humans = administer_drug(humans, y, mda_effectiveness)
 #   else

 #   end

    return humans
end

# function to update the mda information

"""
    update_mda(mda_info, mda_round)

update when the next mda will take place
"""
function update_mda(mda_info, mda_round)

    i = min(mda_round + 1, length(mda_info))
    mda_coverage = mda_info[i].coverage
    min_age_mda =  mda_info[i].min_age
    max_age_mda =  mda_info[i].max_age
    mda_effectiveness =  mda_info[i].effectiveness
    mda_gender = mda_info[i].gender
    if mda_round === length(mda_info)
        next_mda_time = Inf
    else
        next_mda_time = mda_info[i].time
    end
    return mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender
end



#= function to create a set of mda's which will be performed regularly
    first_mda_time specifies when this will first occur in years,
    last_mda_time is the final mda in this block
    regularity is how often to perform the mda in years.
    specify the proportion of pre SAC, SAC and adults at each of these time points
    also specify genders for these differect age groups, along with the effectiveness of mda
=#
"""
    create_mda(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
            last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness)

function to create a set of mda's which will be performed regularly
        first_mda_time specifies when this will first occur in years,
        last_mda_time is the final mda in this block
        regularity is how often to perform the mda in years.
        specify the proportion of pre SAC, SAC and adults at each of these time points
        also specify genders for these differect age groups, along with the effectiveness of mda
"""
function create_mda(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
    last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness)

    mda_info = []
    mda_time = first_mda_time
    while mda_time < last_mda_time
        push!(mda_info, mda_information(pre_SAC_prop, 0, 4, pre_SAC_gender, mda_effectiveness, mda_time))
        push!(mda_info, mda_information(SAC_prop, 4, 16, SAC_gender, mda_effectiveness, mda_time))
        push!(mda_info, mda_information(adult_prop, 16, 110, adult_gender, mda_effectiveness, mda_time))
        mda_time += regularity
    end
    return mda_info
end



# function to add vaccination to population

# function vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
#     vaccine_gender, vaccine_duration, vaccine_round)

#     x = []
#     #= find index of people with correct ages for the mda treatment =#
#     for i in 1 : length(humans)
#         if (humans[i].gender in vaccine_gender) & (min_age_vaccine <= humans[i].age <= max_age_vaccine)
#             push!(x, i)
#         end
#     end


# #=  if this is the first mda round, then treat entirely at random
#     find how many people are eligible for the treatment  =#

#     #if mda_round == 0
#         k = length(x)


# #= randomly permute the indices =#
#         y = shuffle(x)

# #= only take as many as are indicated by the coverage  =#
#         y = y[1:trunc(Int, round(k*vaccine_coverage))]



#         humans = administer_vaccine(humans, y, vaccine_effectiveness, vaccine_duration)
#  #   end

#     #println("output = ", female_worms, male_worms, human_cercariae, eggs)
#     return humans
#     #return x
# end



# #= function to update vaccine information =#
# function update_vaccine(vaccine_info, vaccine_round)

#     i = min(vaccine_round + 1, size(vaccine_info)[1])
#     vaccine_coverage = vaccine_info[i].coverage
#     min_age_vaccine =  vaccine_info[i].min_age
#     max_age_vaccine =  vaccine_info[i].max_age
#     vaccine_duration = vaccine_info[i].duration
#     vaccine_gender = vaccine_info[i].gender
#     if vaccine_round === size(vaccine_info)[1]
#         next_vaccine_time = Inf
#     else
#         next_vaccine_time = vaccine_info[min(vaccine_round + 1, size(vaccine_info)[1])].time
#     end
#     return vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender

# end

"""
    vac_decay!(humans)

decrease vaccination status for each person by 1 each day
"""
function vac_decay!(humans)
    for h in humans
        h.vac_status -= 1
    end
    return humans
end

"""
    kato_katz(eggs, gamma_k)

calculate number of eggs using kato katz method. Gamma_k is a gamma distribution with shape and scale
defined by pars.kato_katz_par
"""
function kato_katz(eggs, gamma_k)
    gamma1 = rand(gamma_k)[1]
    gamma2 = rand(gamma_k)[1]
    pois1 = rand(Poisson(1*gamma1*eggs))[1]
    pois2 = rand(Poisson(1*gamma2*eggs))[1]
    return floor(0.5*(pois1 + pois2))
end



function count_eggs(humans)
    eggs = 0
    for h in humans
        eggs += h.eggs
    end
    return eggs
end






function get_prevalences!(humans, time, pars)

    pop_burden = [0,0,0]
    sac_burden = [0,0,0]
    adult_burden = [0,0,0]
    pop_prev = 0
    sac_prev = 0
    adult_prev = 0
    sac_pop = 0
    adult_pop = 0
    recorded_eggs = Int64[]
    final_ages = Float32[]
    gamma_k  = Gamma(pars.kato_katz_par, 1/kato_katz_par)
    num_humans = length(humans)
    for h in humans
        push!(final_ages, h.age);
        # final_eggs = kato_katz(eggs[i], gamma_k)
        final_eggs = h.eggs * (1-pars.use_kato_katz) + kato_katz(h.eggs, gamma_k) * pars.use_kato_katz
        push!(recorded_eggs, final_eggs)
        if h.age > 5 && h.age < 15
            sac_pop = sac_pop + 1;
        end
        if h.age > 15
            adult_pop = adult_pop + 1;
        end
        if final_eggs > pars.heavy_burden_threshold
            pop_burden[3] = pop_burden[3] + 1
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if 5 < h.age && h.age < 15
                sac_burden[3] = sac_burden[3] + 1
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
            if h.age > 15
                adult_burden[3] = adult_burden[3] + 1
                adult_burden[2] = adult_burden[2] + 1
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        elseif final_eggs > 4
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if h.age > 5 && h.age < 15
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
            if h.age > 15
                adult_burden[2] = adult_burden[2] + 1
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        elseif final_eggs > 0
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if h.age > 5 && h.age < 15
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev +  1
            end
            if h.age > 15
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




#
# repeat_simulations(num_runs, num_sims, N, max_age,
#     initial_worms, worm_stages, max_fecundity, infective_larvae, time_step, contact_rate, r,
#     birth_rate, gamma_pre, female_factor, male_factor)



function save_population_to_file(filename, humans, miracidia, cercariae, pars)
    save(filename, "humans", humans,  "miracidia", miracidia, "cercariae", cercariae, "pars", pars)
end


function load_population_from_file(filename)

    d = load(filename)
    humans = get!(d, "humans",1)
    miracidia = get!(d, "miracidia",1)
    cercariae = get!(d, "cercariae",1)
    pars = get!(d, "pars",1)

    return humans,  miracidia, cercariae, pars
end





#= function to generate a distribution for ages based on a specified demography =#
function generate_age_distribution(pars)

    number_per_age = []
    for i in 1:(length(pars.spec_ages) * pars.ages_per_index)
        index = trunc(Int, (i-1)/pars.ages_per_index) + 1
        push!(number_per_age, pars.spec_ages[index])
    end
    cumsum_spec_ages = cumsum(number_per_age)/sum(number_per_age)
    return cumsum_spec_ages
end


#= function to construct the set of ages, with size N =#
function specified_age_distribution(pars)

    cumsum_spec_ages = generate_age_distribution(pars)
    ages = []
    for i in 1:pars.N
        x = rand()
        k = findall(cumsum_spec_ages .> x)[1]
        push!(ages, k-1)
    end
    return ages
end




function update_env_to_equilibrium(num_time_steps, humans, miracidia, cercariae, pars)



    sim_time = 0
    record_time = 0
    env = Environment()
    record = env.record


    for j in 1:num_time_steps

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        humans =   egg_production!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

#=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)
#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)


    end
    return humans, miracidia, cercariae, record
end





function update_env_to_equilibrium_human_larvae(num_time_steps, humans, miracidia, cercariae, pars)



    sim_time = 0
    record_time = 0
    env = Environment()
    record = env.record


    for j in 1:num_time_steps

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        humans =   egg_production!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

#=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake_with_human_larvae!(humans, cercariae, miracidia, pars)
        humans = human_larvae_maturity(humans, pars)
#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)


    end
    return humans, miracidia, cercariae, record
end



function update_env_to_equilibrium_increasing(num_time_steps, humans, miracidia, cercariae, pars)



    sim_time = 0
    record_time = 0
    env = Environment()
    record = env.record


    for j in 1:num_time_steps

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        humans =   egg_production_increasing!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

#=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)
#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)


    end
    return humans, miracidia, cercariae, record
end






function update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_contact_rate(humans,  pars)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        for h in humans
            h.age += pars.time_step/365

        end

        humans =  egg_production!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

        humans = vac_decay!(humans)

        humans = death_of_human(humans)

        if length(humans) < pars.N
            for k in 1:(pars.N - length(humans))
                humans = birth_of_human(humans, pars)
            end
        end

#=  uptake larvae into humans from the environment  =#
        #=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)


#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =  cercariae_death!(cercariae, pars)

    end
    return humans, miracidia, cercariae, record
end



function update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_contact_rate(humans,  pars)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        for h in humans
            h.age += pars.time_step/365

        end

        humans =  egg_production!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

        humans = vac_decay!(humans)

        humans = death_of_human(humans)

        if length(humans) < pars.N
            for k in 1:(pars.N - length(humans))
                humans = birth_of_human(humans, pars)
            end
        end

#=  uptake larvae into humans from the environment  =#
        #=  uptake larvae into humans from the environment  =#

        humans, cercariae, miracidia = cercariae_uptake_with_human_larvae!(humans, cercariae, miracidia, pars)
        humans = human_larvae_maturity(humans, pars)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =  cercariae_death!(cercariae, pars)

    end
    return humans, miracidia, cercariae, record
end



function update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_contact_rate(humans,  pars)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        for h in humans
            h.age += pars.time_step/365

        end

        humans = egg_production_increasing!(humans, pars)

        humans = worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

        humans = vac_decay!(humans)

        humans = death_of_human(humans)

        if length(humans) < pars.N
            for k in 1:(pars.N - length(humans))
                humans = birth_of_human(humans, pars)
            end
        end

#=  uptake larvae into humans from the environment  =#
        #=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)


#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)

    end
    return humans, miracidia, cercariae, record
end






function update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_contact_rate(humans,  pars)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365



        humans =   egg_production!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

        humans = vac_decay!(humans)

#=  uptake larvae into humans from the environment  =#
        #=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)


#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)

    end
    return humans, miracidia, cercariae, record
end







function update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_contact_rate(humans,  pars)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365



        humans =   egg_production!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

        humans = vac_decay!(humans)

#=  uptake larvae into humans from the environment  =#
        #=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake_with_human_larvae!(humans, cercariae, miracidia, pars)
        humans = human_larvae_maturity(humans, pars)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)

    end
    return humans, miracidia, cercariae, record
end




function update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_contact_rate(humans,  pars)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365



        humans =   egg_production_increasing!(humans, pars)

        humans =  worm_maturity!(humans, pars)

        push!(miracidia, miracidia_production!(humans))

        humans = vac_decay!(humans)

#=  uptake larvae into humans from the environment  =#
        #=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake!(humans, cercariae, miracidia, pars)


#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia =  miracidia_death!(miracidia, pars)
#=  kill cercariae in the environment at specified death rate =#
        cercariae =   cercariae_death!(cercariae, pars)

    end
    return humans, miracidia, cercariae, record
end




# when we run multiple simulations, we store them in an array. This function will store the prevalence and sac prevalence
function collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)
        if run == 1
            for i in 1 : length(record)
                push!(times, record[i].time)
                push!(prev, [record[i].pop_prev])
                push!(sac_prev, [record[i].sac_prev])
                push!(high_burden, [record[i].population_burden[3]])
                push!(high_burden_sac, [record[i].sac_burden[3]])
                push!(adult_prev, [record[i].adult_prev])
                push!(high_adult_burden, [record[i].adult_burden[3]])
            end
        else
            for i in 1 : length(record)
                push!(prev[i], record[i].pop_prev)
                push!(sac_prev[i], record[i].sac_prev)
                push!(high_burden[i], record[i].population_burden[3])
                push!(high_burden_sac[i], record[i].sac_burden[3])
                push!(adult_prev[i], record[i].adult_prev)
                push!(high_adult_burden[i], record[i].adult_burden[3])
            end
        end
    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end



# repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death
function run_repeated_sims_no_population_change(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        humans,  miracidia, cercariae, pars = load_population_from_file(filename)


        humans, miracidia, cercariae, record =
            update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info);



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end




# repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death
function run_repeated_sims_no_population_change_human_larvae(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        humans,  miracidia, cercariae, pars = load_population_from_file(filename)


        humans, miracidia, cercariae, record =
            update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info);



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end




function run_repeated_sims_no_population_change_increasing(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        humans,  miracidia, cercariae, pars = load_population_from_file(filename)


        humans, miracidia, cercariae, record =
            update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info);



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end



# repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death
function run_repeated_sims_no_births_deaths(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        humans,  miracidia, cercariae, pars = load_population_from_file(filename)


        humans, miracidia, cercariae, record =
            update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end



function run_repeated_sims_no_births_deaths_human_larvae(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        humans,  miracidia, cercariae, pars = load_population_from_file(filename)


        humans, miracidia, cercariae, record =
            update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end




# repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death
function run_repeated_sims_no_births_deaths_increasing(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        humans,  miracidia, cercariae, pars = load_population_from_file(filename)


        humans, miracidia, cercariae, record =
            update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,
        high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end



function update_env(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)


    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = pars.record_frequency
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

    for j in 1:num_time_steps

        if sim_time >= update_contact_death_rates
            humans = update_death_rate(humans, death_rate_per_time_step)
            humans = update_contact_rate(humans,  contact_rates_by_age)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences!(humans, sim_time, pars)
            push!(record, a)
            record_time += pars.record_frequency
        end

        sim_time += pars.time_step/365

        for h in humans
            h.age += pars.time_step/365

        end
#=  mature larvae within humans  =#
        humans = human_cercariae_maturity(humans, time_step)

        humans = egg_production(humans, max_fecundity, r, density_dependent_fecundity)

        humans = worm_maturity(humans, worm_stages, average_worm_lifespan, time_step)

        miracidia = miracidia_production(humans, miracidia)

        humans = vac_decay!(humans)

# larvae = larvae_production(d, larvae)
        humans = death_of_human(humans, time_step)

#=  uptake larvae into humans from the environment  =#
        humans, cercariae, miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
#= update information for the next round of mda =#
            mda_round += 1
            mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
                            update_mda(mda_info, mda_round)

        end


#= check if we are at a point in time in which a vaccine is scheduled to take place =#
        if sim_time >= next_vaccine_time

#= perform vaccination =#
            humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                vaccine_gender, vaccine_duration, vaccine_round)
#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                        update_vaccine(vaccine_info, vaccine_round)
        end



#=  kill miracidia in the environment at specified death rate =#
        miracidia = miracidia_death(miracidia, pars)

#=  kill cercariae in the environment at specified death rate =#
        cercariae = cercariae_death(cercariae, env_cercariae_death_rate, time_step)


        l = rand(Binomial(length(humans), birth_rate))[1]
        if l > 0
            for i in 1:l
                humans = birth_of_human(humans, pars)
            end
        end
    end
    return humans, miracidia, cercariae, record
end
