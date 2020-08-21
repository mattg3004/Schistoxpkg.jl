## ##########################
# function to get age dependent death rate.
# The first entry is for under 1's, and the rest are at 5 year intervals
#  At some point this will be changed to be read from a file

# @param {age_death_rate_per_1000} - death rate per 1000 humans per year.
# Hard coded currently
# @param {time step} - how many days we step forward each simulation time step
# @output{death_rate_per_time_step} - age dependent probability of death each time step

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
function generate_ages_and_deaths(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths)
     for i in 1:num_steps
        x = findall(ages - death_ages .> 0)
        k = length(x)
        if k > 0
            for j in k:-1:1
                splice!(ages, x[j])
                splice!(death_ages, x[j])
                push!(ages, 0)
                push!(death_ages, get_death_age(death_prob_by_age, ages_for_deaths))
            end
        end
        ages .+= time_step/365
    end
    return ages, death_ages
end


# function to return the correct death rate per time step based on individuals age

function find_death_rate(age, death_rate_per_time_step)
    if age < 1
        return death_rate_per_time_step[1]
    else
        index = min(2 + trunc(Int, (age-1)/5),length(death_rate_per_time_step))
        return death_rate_per_time_step[index]
    end
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
            contact_rates_by_age = [fill(contact_settings[4], trunc(Int,(max_age+1)))]
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
            contact_rates_by_age = [fill(input_contact_rates[end], trunc(Int,(max_age+1)))]
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
"""
    create_population

    This will create the initial human population with randomly chosen age, and gender.
     Predisposition is taken to be gamma distributed.
     There is also a male and female adjustment to predisposition adjusting for gender specific behaviour
         In addition to this, it will create the initial miracidia environment vector
"""
function create_population(N, max_age, N_communities, community_probs, initial_worms, contact_rates_by_age,
    worm_stages, female_factor, male_factor,
    initial_miracidia, initial_miracidia_days, predis_aggregation, predis_weight,
     time_step,
    mda_adherence, mda_access)

# initialize all the arrays we will keep track of over time
    if length(community_probs) != N_communities
        error("must provide probabilities for membership of each community")
    else
        community_selection = 1
    if N_communities > 1
        community_selection = cumsum(community_probs)/sum(community_probs)
    end

        community =  Int64[]
        female_worms =  Array{Int64}[]
        male_worms = Array{Int64}[]
        human_cercariae = Array{Int64}[]
        eggs =  Int64[]
        vac_status =  Int64[]
        treated =  Int64[]
        vaccinated =  Int64[]
        age_contact_rate =  Float32[]
        death_rate =  Float32[]
        ages =  Float32[]
        gender =  Int64[]
        adherence =  Int64[]
        access =  Int64[]
    #=  initialize the Gamma distribution for predisposition selection  =#
        gamma_pre = Gamma(predis_aggregation, predis_weight/predis_aggregation)

    #=  initialize and fill the environmental variable  =#
        env_miracidia = Int64[]

        for i in 1 : initial_miracidia_days
            push!(env_miracidia, initial_miracidia )
        end

        predisposition = rand(gamma_pre,N)

        for i in 1:N

    #=  begin pushing entries to the data variables we keep track of  =#
            push!(community, findall(community_selection .> rand())[1])
            push!(ages, rand()*max_age)
            push!(gender, rand([0,1]))
            #push!(predisposition, rand(gamma_pre)[1])
            push!(human_cercariae,Int64[])
            push!(eggs,0)
            push!(vac_status, 0)
            push!(treated,0)
            push!(vaccinated, 0)
            if rand() > mda_adherence
                push!(adherence, 0)
            else
                push!(adherence, 1)
            end

            if rand() > mda_access
                push!(access, 0)
            else
                push!(access, 1)
            end
    #=  everyone is initiated with a random number of worms in the first stage  =#
            f_worms = fill(0, worm_stages)
            f_worms[1] = trunc(Int,round(rand()*initial_worms))
            m_worms = fill(0, worm_stages)
            m_worms[1] = trunc(Int,round(rand()*initial_worms))
            push!(female_worms, f_worms)
            push!(male_worms, m_worms)

    #=  age dependent contact rate is found for the given age  =#
            age = (trunc(Int, ages[i]))
            push!(age_contact_rate, contact_rates_by_age[age+1])

    #=  if the person is chosen to be a male of female, then
        adjust their predisposition based on the
        male or female factor, adjusting for
        behavioural differences according to gender  =#
            x = gender[i]
            predisposition[i] = predisposition[i] * (1-x) *female_factor + predisposition[i] * x *female_factor
        end

    #=  return all data that we will use to track the spread of the disease  =#
        return ages , gender, predisposition, community, human_cercariae, eggs, vac_status,
                treated, female_worms, male_worms, vaccinated,
                age_contact_rate,  env_miracidia, adherence, access
    end
end







#= function to update the contact rate of individuals in the population =#
"""
    update_contact_rate(ages, age_contact_rate, contact_rates_by_age)

function to update the contact rate of individuals in the population. This is necessary
    as over time when people age, they will move through different age groups which have
    different contact rates
"""
function update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
    for i in 1:length(ages)
        age = min(length(contact_rates_by_age) - 1, (trunc(Int, ages[i])))
        @inbounds age_contact_rate[i] =  contact_rates_by_age[age+1]
    end
    return age_contact_rate
end








#= function to update all death rates in the population at once =#

function update_death_rate(ages, death_rate, death_rate_per_time_step)

    for i in 1:length(ages)

        if ages[i] < 1
            @inbounds death_rate[i] =  death_rate_per_time_step[1]
        else
            age = trunc(Int, ages[i])
            index = min(2 + trunc(Int, (age-1)/5),length(death_rate_per_time_step))
            @inbounds death_rate[i] = death_rate_per_time_step[index]
        end
    end
    return death_rate
end







# humans uptake larvae based on their predisposition, age_dependent contact rate
# and the number of larvae in the environment. number of larvae taken up is chosen
# from a Poisson distribution
"""
cercariae_uptake(env_miracidia, env_cercariae, time_step, contact_rate,
    community, community_contact_rate, female_worms, male_worms,
    predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop,
    miracidia_maturity_time)

    uptake cer
"""
function cercariae_uptake(env_miracidia, env_cercariae, time_step, contact_rate,
    community, community_contact_rate, female_worms, male_worms,
    predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop,
    miracidia_maturity_time)

#= we want the human population to randomly pick up larvae.
therefore we want to shuffle the population.
the following lines make a random permutation of indices for the population=#
    k = size(female_worms)[1]
    x = randperm(k)
    uptakes = 0
#= assign larvae which have been in the environment for 40 days to become infective.
then delete those larvae from the environmental larvae =#

    if  length(env_miracidia) > (miracidia_maturity_time / time_step)
        env_cercariae += (env_miracidia[1] * human_cercariae_prop)
        splice!(env_miracidia, 1)
    end

#= loop over the population uptaking larvae=#
    for i in 1:k

# set index to the correct value from the random permutation
        @inbounds  j = x[i]
        comm = community[j]
#= if there are still infective larvae in the environment,
we will uptake from choose from Poisson
distribution. otherwise, just uptake 0. =#
    #    if env_cercariae > 0
# calculate the rate of the poisson distribution
            @inbounds  pois_rate  = predisposition[j] * contact_rate * age_contact_rate[j] * community_contact_rate[comm] *
                    env_cercariae * time_step / k

# reduce the rate according to the effectiveness of the vaccine (if any is given)
          # pois_rate = pois_rate * (1 - (vac_status[j] > 0) * vaccine_effectiveness);

# choose from the Poisson distribution
            uptake = rand(Poisson(pois_rate))
    #    else
    #        uptake = 0
    #    end
        if uptake >0
# push the number of uptaken larvae into the human larvae array
            n = rand(Binomial(uptake, 0.5))

            female_worms[j][1] += n
            male_worms[j][1] += uptake - n
        end

        # push!(human_cercariae[j], uptake)

# # reduce the infective larvae by the number of larvae uptaken
        env_cercariae -= uptake

        end

# return the infective, human and environmental larvae arrays
    return env_cercariae, env_miracidia, female_worms, male_worms

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
function cercariae_death(env_cercariae, env_cercariae_survival_prop, time_step)
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


function worm_maturity(female_worms, male_worms, worm_stages,
    average_worm_lifespan, time_step)



# probability of aging out of category/ dying
    p = time_step * worm_stages / ( 365 * average_worm_lifespan)

# loop over the worms
    for i in 1:size(female_worms)[1]

# kill appropriate number of worms in the final stage
        n = female_worms[i][worm_stages]
        # if n > 0
            # dis = Binomial(n, 1-p)
            female_worms[i][worm_stages] = rand(Binomial(n, 1-p))
        # end


        n = male_worms[i][worm_stages]
        # if n > 0
            # dis = Binomial(n, 1-p)
            male_worms[i][worm_stages] = rand(Binomial(n, 1-p))
        # end

#=
     for aging worms, we do this in reverse order, which ensures the
     correct order of aging is respected
=#

        for j in (worm_stages-1):-1:1
#=   choose the number of male and female worms to age from one stage to the next   =#
            aging_females = rand(Binomial(female_worms[i][j], p))
            aging_males = rand(Binomial(male_worms[i][j], p))

#=   add and subtract the number of worms from the appropriate categories   =#
            female_worms[i][j+1] += aging_females
            female_worms[i][j] -= aging_females
            male_worms[i][j+1] += aging_males
            male_worms[i][j] -= aging_males
        end
    end

#=  return the female and male worm arrays  =#
    return female_worms, male_worms

end



# function to calculate the number of male-female pairs of worms in each human
# this is just the minimum of total female worms and the total male worms


function calculate_worm_pairs(female_worms, male_worms)
    worm_pairs = Int64[]
    for i in 1:length(female_worms)
        push!(worm_pairs, min(sum(female_worms[i]), sum(male_worms[i])))
    end
    return(worm_pairs)
    # return min.(sum.(female_worms), sum.(male_worms))
end





# function to calculate the number of male and female worms in each human

function calculate_total_worms(female_worms, male_worms)
    return sum.(female_worms), sum.(male_worms)
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
function egg_production(eggs, max_fecundity, r, worm_pairs,
                        density_dependent_fecundity, time_step)


# loop over individuals
    for i in 1 : length(eggs)

#= if we have a positive number of worms, then make calculation,
otherwise the number of eggs is trivially 0 =#
        #@inbounds if worm_pairs[i] > 0

# calculate the mean number of eggs we would expect
                # @inbounds    mean_eggs =  max_fecundity * worm_pairs[i] *
                #     exp(- density_dependent_fecundity *
                #     (total_female_worms[i] + total_male_worms[i]))
wp = worm_pairs[i]
wp = max(wp,0.000000001)
                    @inbounds    mean_eggs =  max_fecundity * wp *
                        exp(- density_dependent_fecundity *
                        (wp ))

# calculate the number of successes
                @inbounds      NB_r = r * wp

# calculate the probability of a success
                p = NB_r/(NB_r+mean_eggs)

# choose from NB

                eggs_num = rand(NegativeBinomial(NB_r,p))[1]

        #    else
        #        eggs_num = 0
        #    end

# put this selected number of eggs into the eggs array
            @inbounds eggs[i] = eggs_num
        end

# return the eggs array
    return eggs
end


# hatch the eggs in the humans into the environment

function miracidia_production(eggs, env_miracidia, time_step, age_contact_rate, community_contact_rate, community)
#= as we can step forward an arbitrary number of days at a time, we multiply the number of miracidia by the
    length of the forward step, assuming that each of the last given number of days were equivalent to each other
=#
    max_contact_rate = maximum(age_contact_rate)
    xx = age_contact_rate ./ max_contact_rate
    released_eggs = 0.0
    for i in 1:length(eggs)
        released_eggs +=  xx[i] * eggs[i] * community_contact_rate[community[i]]/maximum(community_contact_rate)
    end
    push!(env_miracidia,  round(sum(released_eggs)))
    return env_miracidia
end





# function to kill humans at age dependent rate

function death_of_human(ages, death_ages, gender, predisposition,  human_cercariae, eggs,
                            vac_status, treated, female_worms, male_worms,
                            vaccinated, age_contact_rate, time_step,
                            adherence, access)

#= loop through the population, and based on the death rate,
delete individuals from the population  =#

        time_to_death = death_ages .- ages

#= if there is anyone who is older than their death age, then find them =#
        x = findall(time_to_death .< 0)
        k = length(x)
        if k > 0
            for j in k:-1:1
                i = x[j]
                splice!(ages, i)
                splice!(death_ages, i)
                splice!(gender, i)
                splice!(predisposition, i)
                splice!(human_cercariae, i)
                splice!(eggs, i)
                splice!(vac_status, i)
                splice!(treated, i)
                splice!(female_worms, i)
                splice!(male_worms, i)
                splice!(vaccinated, i)
                splice!(age_contact_rate, i)
                splice!(adherence, i)
                splice!(access, i)
            end
        end


#=  return the arrays  =#
    return ages, death_ages, gender, predisposition,  human_cercariae, eggs,
                                vac_status, treated, female_worms, male_worms,
                                vaccinated, age_contact_rate,
                                adherence, access
end







# function to add a person to the data if a birth occurs

function birth_of_human(ages, death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
                        treated, female_worms, male_worms,vaccinated, age_contact_rate,
                        female_factor, male_factor, contact_rates_by_age,
                        worm_stages, predis_aggregation, predis_weight,
                        adherence, death_prob_by_age, ages_for_deaths,community_probs,
                        mda_adherence, access, mda_access)


#  load the gamma distribution for the predispostion distribution
    gamma_pre = Gamma(predis_aggregation, predis_weight/predis_aggregation)


# fill female and male worms array
    f_worms = fill(0, worm_stages)
    m_worms = fill(0, worm_stages)
    community_selection = 1
    if length(community_probs) > 1
        community_selection = cumsum(community_probs)/sum(community_probs)
    end
    push!(community, findall(community_selection .> rand())[1])

# push into arrays
    push!(female_worms, f_worms)
    push!(male_worms, m_worms)
    push!(ages, 0)
    push!(gender, rand([0,1]))
    push!(predisposition, rand(gamma_pre)[1])
    push!(human_cercariae, Int64[])
    push!(eggs,0)
    push!(vac_status, 0)
    push!(treated, 0)
    push!(vaccinated, 0)
    predisp = rand(gamma_pre)[1]

    push!(death_ages, get_death_age(death_prob_by_age, ages_for_deaths))
    push!(age_contact_rate, contact_rates_by_age[1])
    if rand() > mda_adherence
        push!(adherence, 0)
    else
        push!(adherence, 1)
    end

    if rand() > mda_access
        push!(access, 0)
    else
        push!(access, 1)
    end
# adjust predisposition based on gender specific behaviour parameters
    if gender[end] === 0
        predisposition[end] = predisp * female_factor
    else
        predisposition[end] = predisp * male_factor
    end

#  return arrays
    return ages, death_ages, gender, predisposition, community,
    human_cercariae, eggs, vac_status,
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, adherence, access
end










# mature the larvae within humans into worms after 35 days


function human_cercariae_maturity(human_cercariae, female_worms, male_worms, time_step)

#=  loop over humans  =#
    for i in 1:size(human_cercariae)[1]

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_cercariae  =#
        if length(human_cercariae[i]) > round(35/time_step, digits = 0)
            females = rand(Binomial(human_cercariae[i][1], 0.5))[1]
            @inbounds female_worms[i][1] += females
            @inbounds male_worms[i][1] += human_cercariae[i][1] - females
            @inbounds splice!(human_cercariae[i],1)
        end
    end

#= return arrays  =#
    return human_cercariae, female_worms, male_worms
end








# function to administer drug to a specific variable (e.g. female_worms or eggs).
# input the variable, the indices to apply to and the effectiveness of treatment

function administer_drug(d, indices, drug_effectiveness, adherence)
    if drug_effectiveness === 1
        @inbounds d[indices] .*= (1 .- adherence[indices])

    elseif drug_effectiveness === 0

    else
        for i in 1:length(indices)
            @inbounds index = indices[i]
            d[index] = trunc.(Int, d[index])
            @inbounds d[index] = rand.(Binomial.(d[index], 1 - (drug_effectiveness * adherence[index])))
        end
    end
    return d
end







# function for mass drug administration
# currently there is no correlation between individuals chosen each time

function mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender,
    ages, female_worms, male_worms, human_cercariae, eggs,
    treated, mda_round, gender, adherence, access)

#= find index of people with correct ages for the mda treatment =#
    in_gender = in(mda_gender)
    x = findall(((min_age_mda .<= ages .<= max_age_mda) .& in_gender.(gender) .& (access .== 1)))

#=  if this is the first mda round, then treat entirely at random
    find how many people are eligible for the treatment  =#

    #if mda_round == 0
        k = length(x)


#= randomly permute the indices =#
        y = shuffle(x)

#= only take as many as are indicated by the coverage  =#
        y = y[1:trunc(Int, round(k*mda_coverage))]

# update female and male worms, human cercariae and eggs
        female_worms = administer_drug(female_worms, y, mda_effectiveness, adherence)
        male_worms = administer_drug(male_worms, y, mda_effectiveness, adherence)
    #    human_cercariae = administer_drug(human_cercariae, y, mda_effectiveness, adherence)
        eggs = administer_drug(eggs, y, 1, adherence)
 #   else

 #   end

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
    if mda_round >= size(mda_info)[1]
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
function create_mda(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
    last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness)

    mda_info = []
    mda_time = first_mda_time
    while mda_time <= last_mda_time
        push!(mda_info, mda_information(pre_SAC_prop, 0, 5, pre_SAC_gender, mda_effectiveness, mda_time))
        push!(mda_info, mda_information(SAC_prop, 5, 15, SAC_gender, mda_effectiveness, mda_time))
        push!(mda_info, mda_information(adult_prop, 15, 110, adult_gender, mda_effectiveness, mda_time))
        mda_time += regularity
    end
    return mda_info
end



# function to add vaccination to population

function vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
    vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
    treated, vaccine_duration, vac_status, vaccine_round, gender, access)

#= find index of people with correct ages for the mda treatment =#
    # x = findall( min_age_vaccine .<= ages .<= max_age_vaccine)

    in_gender = in(vaccine_gender)
    x = findall(((min_age_vaccine .<= ages .<= max_age_vaccine) .& in_gender.(gender).& (access .== 1)))
#=  if this is the first mda round, then treat entirely at random
    find how many people are eligible for the treatment  =#

    #if mda_round == 0
        k = length(x)
        if k > 0
            y = shuffle(x)

    #= only take as many as are indicated by the coverage  =#
            y = y[1:trunc(Int, round(k*vaccine_coverage))]

    #= since vaccination is administered by a health care worker, there is no issue of adherence,
    therefore make sure that all people adhere to the treatment for vaccination =#
            ad = ones(length(female_worms))
    # update female and male worms, human cercariae and eggs post vaccine
            female_worms = administer_drug(female_worms, y, vaccine_effectiveness, ad)
            male_worms = administer_drug(male_worms, y, vaccine_effectiveness, ad)
            human_cercariae = administer_drug(human_cercariae, y, vaccine_effectiveness, ad)
            eggs = administer_drug(eggs, y, 1, ad)
            vac_status[y] .= vaccine_duration

        end


    return female_worms, male_worms, human_cercariae, eggs, vac_status
end






#= function to update vaccine information =#
function update_vaccine(vaccine_info, vaccine_round)

    i = min(vaccine_round + 1, size(vaccine_info)[1])
    vaccine_coverage = vaccine_info[i].coverage
    min_age_vaccine =  vaccine_info[i].min_age
    max_age_vaccine =  vaccine_info[i].max_age
    vaccine_duration = vaccine_info[i].duration
    vaccine_gender = vaccine_info[i].gender
    if vaccine_round >= size(vaccine_info)[1]
        next_vaccine_time = Inf
    else
        next_vaccine_time = vaccine_info[min(vaccine_round + 1, size(vaccine_info)[1])].time
    end
    return vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender

end


# function to update the variable arrays a given number of times (num_time_steps)

function update_env(num_time_steps, ages, death_ages, community, community_contact_rate, community_probs,
    human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,predis_weight,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
    vaccinated, age_contact_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
    record_frequency, human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold,
    kato_katz_par, use_kato_katz)

    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = record_frequency
    record = out[]
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
            age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences(ages, eggs, time, heavy_burden_threshold, kato_katz_par, use_kato_katz)
            push!(record, a)
            record_time += record_frequency
        end

        sim_time += time_step/365
#=  add to age variables  =#
        ages = ages .+ time_step/365

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

                               density_dependent_fecundity, time_step)

#=  mature worms in each human  =#
        female_worms, male_worms = worm_maturity(female_worms, male_worms,
                                                 worm_stages, average_worm_lifespan,
                                                 time_step)

 #=  reduce the vaccination status by the time step  =#
        vac_status = vac_status .- time_step/365

#=  hacth the human eggs into the environment  =#
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step, age_contact_rate, community_contact_rate, community)

#=  update population due to death  =#
        ages, death_ages, gender, predisposition,  human_cercariae, eggs,
                                    vac_status, treated, female_worms, male_worms,
                                    vaccinated, age_contact_rate,
                                    adherence, access =
            death_of_human(ages, death_ages, gender, predisposition,  human_cercariae, eggs,
                                        vac_status, treated, female_worms, male_worms,
                                        vaccinated, age_contact_rate,  time_step,
                                        adherence, access)


#=  uptake larvae into humans from the environment  =#
        env_cercariae, env_miracidia, female_worms, male_worms =
        cercariae_uptake(env_miracidia, env_cercariae, time_step, contact_rate,
            community, community_contact_rate, female_worms, male_worms,
            predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop,
            miracidia_maturity_time)


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
                    treated, vaccine_duration, vac_status, vaccine_round, gender, access)

#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                update_vaccine(vaccine_info, vaccine_round)
        end

#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_survival_prop)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_survival_prop, time_step)

 #=  choose from binomial distribution for the number of births in the population  =#
        l = rand(Binomial(size(ages)[1], birth_rate))[1]

 #=  loop over the number of births that there are  =#
        if l > 0
            for i in 1:l

 #=  update the population due to births  =#
                 ages, death_ages, gender, predisposition, community,
                 human_cercariae, eggs, vac_status,
                 treated, female_worms, male_worms,
                 vaccinated, age_contact_rate, adherence, access =
                     birth_of_human(ages, death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
                                             treated, female_worms, male_worms,vaccinated, age_contact_rate,
                                             female_factor, male_factor, contact_rates_by_age,
                                             worm_stages, predis_aggregation, predis_weight,
                                             adherence, death_prob_by_age, ages_for_deaths,community_probs,
                                             mda_adherence, access, mda_access)
            end
        end
    end

#=  return the arrays  =#
    return ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate,
    env_miracidia, env_cercariae, adherence,access,
    record
end





# function to perform kato katz procedure to count eggs in stool
function kato_katz(eggs, gamma_k)
    gamma1 = rand(gamma_k)[1]
    gamma2 = rand(gamma_k)[1]
    pois1 = rand(Poisson(1*gamma1*eggs))[1]
    pois2 = rand(Poisson(1*gamma2*eggs))[1]
    return floor(0.5*(pois1 + pois2))

end



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


function get_prevalences(ages, eggs, time, heavy_burden_threshold, kato_katz_par, use_kato_katz)

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

    num_humans = size(ages)[1]
    gamma_k = Gamma(kato_katz_par, 1 / kato_katz_par)
    for i in 1:num_humans
        push!(final_ages, ages[i]);
        final_eggs = use_kato_katz * kato_katz(eggs[i], gamma_k) + (1-use_kato_katz)* eggs[i]
        # final_eggs = eggs[i]
        push!(recorded_eggs, final_eggs)
        if ages[i] > 5 && ages[i] < 15
            sac_pop = sac_pop + 1;
        end
        if ages[i] > 15
            adult_pop = adult_pop + 1;
        end
        if final_eggs > heavy_burden_threshold
            pop_burden[3] = pop_burden[3] + 1
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if 5 < ages[i] && ages[i] < 15
                sac_burden[3] = sac_burden[3] + 1
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
            if ages[i] > 15
                adult_burden[3] = adult_burden[3] + 1
                adult_burden[2] = adult_burden[2] + 1
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        elseif final_eggs > 4
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if ages[i] > 5 && ages[i] < 15
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
            if ages[i] > 15
                adult_burden[2] = adult_burden[2] + 1
                adult_burden[1] = adult_burden[1] + 1
                adult_prev = adult_prev + 1
            end
        elseif final_eggs > 0
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if ages[i] > 5 && ages[i] < 15
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev +  1
            end
            if ages[i] > 15
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


function load_population_from_file(filename, chosen_pop_size, full_size)
    d = load(filename);
    ages = get!(d, "ages",1)
    gender = get!(d, "gender",1)
    predisposition = get!(d, "predisposition",1)
    human_cercariae = get!(d, "human_cercariae",1)
    eggs = get!(d, "eggs",1)
    vac_status = get!(d, "vac_status",1)
    treated = get!(d, "treated",1)
    female_worms = get!(d, "female_worms",1)
    male_worms = get!(d, "male_worms",1)
    vaccinated = get!(d, "vaccinated",1)
    age_contact_rate = get!(d, "age_contact_rate",1)
    death_ages = get!(d, "death_ages",1)
    env_miracidia = get!(d, "env_miracidia",1)
    env_cercariae = get!(d, "env_cercariae",1)
    adherence = get!(d, "adherence",1)
    access = get!(d, "access",1)
    community = get!(d, "community",1)

    if full_size != true
        x = randperm(length(ages))[1:choose_pop_size];

        ages = ages[x]
        gender = gender[x]
        predisposition = predisposition[x]
        human_cercariae = human_cercariae[x]
        eggs = eggs[x]
        vac_status = vac_status[x]
        treated = treated[x]
        female_worms = female_worms[x]
        male_worms = male_worms[x]
        vaccinated = vaccinated[x]
        age_contact_rate = age_contact_rate[x]
        death_ages = death_ages[x]
        adherence = adherence[x]
        access = access[x]
        community = community[x]
    end
     return ages, death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status, treated,
        female_worms, male_worms, vaccinated, age_contact_rate,
        env_miracidia, env_cercariae, adherence, access
end



function save_population_to_file(filename, ages, gender, predisposition, community, human_cercariae, eggs, vac_status, treated,
        female_worms, male_worms, vaccinated, age_contact_rate, death_ages, env_miracidia, env_cercariae, adherence, access)

    save(filename, "ages", ages ,  "gender", gender,"predisposition",   predisposition,
        "human_cercariae", human_cercariae, "eggs", eggs, "community", community,
        "vac_status", vac_status,"treated", treated, "female_worms",  female_worms, "male_worms", male_worms,
        "vaccinated", vaccinated,  "age_contact_rate", age_contact_rate, "death_ages", death_ages,
        "env_miracidia",env_miracidia, "env_cercariae", env_cercariae, "adherence", adherence, "access", access)
end



function run_simulation_from_loaded_population(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    community, community_contact_rate, community_probs,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info,
    record_frequency)


    ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate,
    env_miracidia, env_cercariae, adherence,access,
    record =
    update_env(num_time_steps, ages, death_ages, community, community_contact_rate, community_probs,
        human_cercariae, female_worms, male_worms,
        time_step, average_worm_lifespan,
        eggs, max_fecundity, r, worm_stages,
        vac_status, gender, predis_aggregation,predis_weight,
        predisposition, treated, vaccine_effectiveness,
        density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
        vaccinated, age_contact_rate, env_miracidia,
        env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
        female_factor, male_factor, contact_rates_by_age,
        birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
        record_frequency, human_cercariae_prop,miracidia_maturity_time)

    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms, vaccinated,
    age_contact_rate, death_rate, env_miracidia, env_cercariae, record, adherence, access

end




function find_mean_worms_by_age(ages, female_worms, male_worms, sim_number, mean_worms_by_age)
    rounded_ages = trunc.(Int,round.(ages, digits = 0))
    if sim_number == 1
        mean_worms_by_age =[]
        for i in 0:maximum(rounded_ages)
            x = findall(rounded_ages .== i)
            if length(x) > 0
                age_worms = sum(sum.(female_worms[x]) + sum.(male_worms[x]))
                push!(mean_worms_by_age, [i, mean(age_worms), length(x)])
            end
        end
    else
        for i in 0:maximum(rounded_ages)
            index = trunc(Int, i+1)
            if index <= length(mean_worms_by_age)
                x = findall(rounded_ages .== i)
                if length(x) > 0
                    age_worms = sum(sum.(female_worms[x]) + sum.(male_worms[x]))
                    mean_worms_by_age[index][2] += age_worms
                    mean_worms_by_age[index][3] += length(x)
                end
            end
        end
    end
    return mean_worms_by_age
end


function run_simulation(N, max_age, initial_worms, time_step, worm_stages, female_factor, male_factor,
    initial_miracidia, initial_miracidia_days,env_cercariae, contact_rate, max_fecundity,
    density_dependent_fecundity, r, num_time_steps, birth_rate, average_worm_lifespan, predis_aggregation,
    env_cercariae_survival_prop, env_miracidia_survival_prop, mda_coverage, mda_round, vaccine_effectiveness,
    mda_info, vaccine_info, record_frequency, mda_adherence, scenario,human_cercariae_prop)


    age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                               0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                               21.83, 29.98, 36.98]


    contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [], [])
    death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)
        death_rate_per_time_step =     death_rate_per_time_step.*100


    ages , gender, predisposition,  human_cercariae, eggs, vac_status,
            treated, female_worms, male_worms, vaccinated,
            age_contact_rate, env_miracidia, adherence, access =
        create_population(N, max_age, N_communities, community_probs, initial_worms, contact_rates_by_age,
            worm_stages, female_factor, male_factor,
            initial_miracidia, initial_miracidia_days, predis_aggregation, predis_weight,
             time_step,
            mda_adherence, mda_access)



            ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
            vac_status, treated, female_worms, male_worms,
            vaccinated, age_contact_rate,
            env_miracidia, env_cercariae, adherence,access,
            record =
        update_env(num_time_steps, ages, death_ages, community, community_contact_rate, community_probs,
            human_cercariae, female_worms, male_worms,
            time_step, average_worm_lifespan,
            eggs, max_fecundity, r, worm_stages,
            vac_status, gender, predis_aggregation,predis_weight,
            predisposition, treated, vaccine_effectiveness,
            density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
            vaccinated, age_contact_rate, env_miracidia,
            env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
            female_factor, male_factor, contact_rates_by_age,
            birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
            record_frequency, human_cercariae_prop,miracidia_maturity_time)

    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms, vaccinated,
    age_contact_rate, death_rate, env_miracidia, env_cercariae,
    adherence, record
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################




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



#= same function to create the population as create_population
    but with a pre-specified distribution of ages
=#
"""
    create_population_specified_ages

This will create the initial human population with an age distribution
specified by the spec_ages variable
Predisposition is taken to be gamma distributed. There is also a male and female
adjustment to predisposition adjusting for gender specific behaviour
In addition to this, it will create the initial miracidia environment vector
"""
function create_population_specified_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
        worm_stages, female_factor, male_factor,initial_miracidia,
        initial_miracidia_days, predis_aggregation, predis_weight,
        time_step,
        spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
        mda_adherence, mda_access)

    if length(community_probs) != N_communities
        error("must provide probabilities for membership of each community")
    else
        community_selection = 1
    if N_communities > 1
        community_selection = cumsum(community_probs)/sum(community_probs)
    end
        community =  Int64[]
    # initialize all the arrays we will keep track of over time

        female_worms = Array{Int64}[]
        male_worms = Array{Int64}[]
        human_cercariae = Array{Int64}[]
        eggs = Int64[]
        vac_status = Int64[]
        treated = Int64[]
        vaccinated = Int64[]
        age_contact_rate = Float32[]
        death_rate = Float32[]
        gender = Int64[]
        adherence = Int64[]
        access = Int64[]
        death_ages = Float32[]
    #=  initialize the Gamma distribution for predisposition selection  =#
        gamma_pre = Gamma(predis_aggregation, predis_weight/predis_aggregation)

    #=  initialize and fill the environmental variable  =#
        env_miracidia = Int64[]

        for i in 1 : initial_miracidia_days
            push!(env_miracidia, initial_miracidia )
        end

        predisposition = rand(gamma_pre,N)
        ages = specified_age_distribution(N, spec_ages, ages_per_index)

        for i in 1:N

    #=  begin pushing entries to the data variables we keep track of  =#
            push!(community, findall(community_selection .> rand())[1])
            push!(death_ages, get_death_age(death_prob_by_age, ages_for_deaths))
            push!(gender, rand([0,1]))
            push!(human_cercariae,Int64[])
            push!(eggs,0)
            push!(vac_status, 0)
            push!(treated,0)
            push!(vaccinated, 0)
            if rand() > mda_adherence
                push!(adherence, 0)
            else
                push!(adherence, 1)
            end
            if rand() > mda_access
                push!(access, 0)
            else
                push!(access, 1)
            end
    #=  everyone is initiated with a random number of worms in the first stage  =#
            f_worms = fill(0, worm_stages)
            f_worms[1] = round(rand()*initial_worms)
            m_worms = fill(0, worm_stages)
            m_worms[1] = round(rand()*initial_worms)
            push!(female_worms, f_worms)
            push!(male_worms, m_worms)

    #=  age dependent contact rate is found for the given age  =#
            age = (trunc(Int, ages[i]))
            push!(age_contact_rate, contact_rates_by_age[age+1])
            ages[i] += rand()

    #=  if the person is chosen to be a male of female, then
        adjust their predisposition based on the
        male or female factor, adjusting for
        behavioural differences according to gender  =#
            x = gender[i]
            predisposition[i] = predisposition[i] * (1-x) *female_factor + predisposition[i] * x *female_factor
        end

    #=  return all data that we will use to track the spread of the disease  =#
        return ages , death_ages, gender, predisposition, community,
            human_cercariae, eggs, vac_status,
            treated, female_worms, male_worms, age_contact_rate,
            vaccinated, env_miracidia, adherence, access
    end
end



#= this function will run the population for a specified number of time steps
with no births, deaths or aging in the population, hence running it to equilibrium
for the specified size of population =#

function update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    community, community_contact_rate,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,vaccinated, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate,human_cercariae_prop,
    miracidia_maturity_time, heavy_burden_threshold, kato_katz_par, use_kato_katz)


    sim_time = 0
    record_time = record_frequency
    record = out[]


#=  loop for number of sims  =#

    for j in 1:num_time_steps


        if sim_time >= record_time
            a = get_prevalences(ages, eggs, time, heavy_burden_threshold, kato_katz_par, use_kato_katz)
            push!(record, a)
            record_time += record_frequency
        end

        sim_time += time_step/365


#=  mature larvae within humans  =#
# print("human_cercariae_maturity")
#    @time human_cercariae, female_worms, male_worms =
#             human_cercariae_maturity(human_cercariae, female_worms, male_worms, time_step)


#=  calculate the number of worm pairs in each human  =#
    worm_pairs = calculate_worm_pairs(female_worms, male_worms)

#=  calculate the total number of worms in each human  =#
# print("calculate total worms")
#        total_female_worms, total_male_worms =
#     @time    calculate_total_worms(female_worms, male_worms)
    eggs =  egg_production(eggs, max_fecundity, r, worm_pairs,
                               density_dependent_fecundity, time_step)

#=  mature worms in each human  =#
    female_worms, male_worms = worm_maturity(female_worms, male_worms,
                                                 worm_stages, average_worm_lifespan,
                                                 time_step)



#=  hacth the human eggs into the environment  =#
    env_miracidia =  miracidia_production(eggs, env_miracidia, time_step, age_contact_rate, community_contact_rate, community)

#=  uptake larvae into humans from the environment  =#
    env_cercariae, env_miracidia, female_worms, male_worms =
        cercariae_uptake(env_miracidia, env_cercariae, time_step, contact_rate,
            community, community_contact_rate, female_worms, male_worms,
            predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop,
            miracidia_maturity_time)



#=  kill miracidia in the environment at specified death rate =#
       env_miracidia = miracidia_death(env_miracidia, env_miracidia_survival_prop)

#=  kill cercariae in the environment at specified death rate =#
    env_cercariae = cercariae_death(env_cercariae, env_cercariae_survival_prop, time_step)

    end

#=  return the arrays  =#
    return ages, gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, env_miracidia, env_cercariae, record
end




# allow mda and vaccine when updating the model, but whenever there is a death, we instantly add another individual
function update_env_keep_population_same(num_time_steps, ages, death_ages,community, community_contact_rate, community_probs,
    human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,predis_weight,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
    vaccinated, age_contact_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
    record_frequency, human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold,
    kato_katz_par, use_kato_katz)

    start_pop = length(ages)
    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = 0
    record = out[]
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
            age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
            update_contact_death_rates += 1/5
        end

        if sim_time >= record_time
            a = get_prevalences(ages, eggs, time, heavy_burden_threshold, kato_katz_par, use_kato_katz)
            push!(record, a)
            record_time += record_frequency
        end

        sim_time += time_step/365

        ages .+= time_step/365
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

                               density_dependent_fecundity, time_step)

#=  mature worms in each human  =#
        female_worms, male_worms = worm_maturity(female_worms, male_worms,
                                                 worm_stages, average_worm_lifespan,
                                                 time_step)

 #=  reduce the vaccination status by the time step  =#
        vac_status = vac_status .- time_step/365

#=  hacth the human eggs into the environment  =#
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step, age_contact_rate, community_contact_rate, community)

    #=  update population due to death  =#
    ages, death_ages, gender, predisposition,  human_cercariae, eggs,
                                vac_status, treated, female_worms, male_worms,
                                vaccinated, age_contact_rate,
                                adherence, access =
        death_of_human(ages, death_ages, gender, predisposition,  human_cercariae, eggs,
                                    vac_status, treated, female_worms, male_worms,
                                    vaccinated, age_contact_rate,  time_step,
                                    adherence, access)

        if length(ages) < start_pop
            num_births = start_pop - length(ages)
            for i in 1:num_births
                ages, death_ages, gender, predisposition, community,
                human_cercariae, eggs, vac_status,
                treated, female_worms, male_worms,
                vaccinated, age_contact_rate, adherence, access =
                    birth_of_human(ages, death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
                                            treated, female_worms, male_worms,vaccinated, age_contact_rate,
                                            female_factor, male_factor, contact_rates_by_age,
                                            worm_stages, predis_aggregation, predis_weight,
                                            adherence, death_prob_by_age, ages_for_deaths,community_probs,
                                            mda_adherence, access, mda_access)
            end
        end
#=  uptake larvae into humans from the environment  =#
        env_cercariae, env_miracidia, female_worms, male_worms =
        cercariae_uptake(env_miracidia, env_cercariae, time_step, contact_rate,
            community, community_contact_rate, female_worms, male_worms,
            predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop,
            miracidia_maturity_time)


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
                    treated, vaccine_duration, vac_status, vaccine_round, gender, access)

#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                update_vaccine(vaccine_info, vaccine_round)
        end

#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_survival_prop)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_survival_prop, time_step)
        # println("time = ", sim_time)
        # println( "female_worms =", sum(sum(female_worms)))
        # println( "male_worms =", sum(sum(male_worms)))
        # println( "age_contact_rate =", sum(age_contact_rate))
        # println("cerc = ", env_cercariae)
        # println("mira = ", env_miracidia)
    end

#=  return the arrays  =#
    return ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate,
    env_miracidia, env_cercariae, adherence,access,
    record
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
function run_repeated_sims_no_population_change(num_repeats, num_time_steps,
    time_step, average_worm_lifespan, community_contact_rate, community_probs,
    max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
    density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
    record_frequency, filename, human_cercariae_prop, miracidia_maturity_time,  heavy_burden_threshold)

    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac =[]
    adult_prev = []
    high_adult_burden = []


    for run in 1:num_repeats


        ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
        human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
         env_cercariae_equ, adherence_equ, access_equ =
        load_population_from_file(filename, N, true)

        ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate,
        env_miracidia, env_cercariae, adherence,access,
        record =
         update_env_keep_population_same(num_time_steps, copy(ages_equ), copy(death_ages_equ),
         copy(community_equ), community_contact_rate, community_probs,
         copy(human_cercariae_equ), copy(female_worms_equ), copy(male_worms_equ),
                    time_step, average_worm_lifespan,
                    copy(eggs_equ), max_fecundity, r, worm_stages,
                    copy(vac_status_equ), copy(gender_equ), predis_aggregation, predis_weight,
                    copy(predisposition_equ), copy(treated_equ), vaccine_effectiveness,
                    density_dependent_fecundity,death_prob_by_age, ages_for_deaths,
                    copy(vaccinated_equ) , copy(age_contact_rate_equ), copy(env_miracidia_equ) ,
                    copy(env_cercariae_equ) , contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                    female_factor, male_factor, contact_rates_by_age,
                    birth_rate, mda_info, vaccine_info, copy(adherence_equ), mda_adherence,
                    copy(access_equ), mda_access,
                    record_frequency,human_cercariae_prop,miracidia_maturity_time, heavy_burden_threshold);



        times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden =
            collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)

    end

    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end




# repeat simulations where we allow mdas and vaccination, and allow random births and deaths
function run_repeated_sims_random_births_deaths(num_repeats, num_time_steps,
                time_step, average_worm_lifespan,
                max_fecundity, r, worm_stages,
                predis_aggregation,predis_weight,
                density_dependent_fecundity,
                age_contact_rate_equ, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                female_factor, male_factor, contact_rates_by_age,death_prob_by_age, ages_for_deaths,
                birth_rate, mda_info, vaccine_info,  mda_adherence, mda_access,
                record_frequency, times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, filename,
                human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold)


                times = []
                prev = []
                sac_prev = []
                high_burden = []
                high_burden_sac =[]
                adult_prev = []
                high_adult_burden = []

    for run in 1:num_repeats

        ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
         env_cercariae_equ, adherence_equ, access_equ =
        load_population_from_file(filename, N, true)




        ages, death_ages, gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate,
        env_miracidia, env_cercariae, adherence,access,
        record =
         update_env(num_time_steps, copy(ages_equ), copy(death_ages_equ),copy(human_cercariae_equ), copy(female_worms_equ), copy(male_worms_equ),
                    time_step, average_worm_lifespan,
                    copy(eggs_equ), max_fecundity, r, worm_stages,
                    copy(vac_status_equ), copy(gender_equ), predis_aggregation,predis_weight,
                    copy(predisposition_equ), copy(treated_equ), vaccine_effectiveness,
                    density_dependent_fecundity,death_prob_by_age, ages_for_deaths,
                    copy(vaccinated_equ) , copy(age_contact_rate_equ),  copy(env_miracidia_equ) ,
                    copy(env_cercariae_equ) , contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                    female_factor, male_factor, contact_rates_by_age,
                    birth_rate, mda_info, vaccine_info, copy(adherence_equ), mda_adherence,
                    copy(access_equ), mda_access,
                    record_frequency,human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold);

                    times, prev, sac_prev, high_burden, high_burden_sac, adult_prev = collect_prevs(times, prev, sac_prev, high_burden,
                    high_burden_sac, adult_prev, high_adult_burden, record, run)

    end
    return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden
end



function plot_sac_burden_and_sac_high_burden(r)
    # times = Float32[]
    # prev = Float32[]
    # sac_prev = Float32[]
    # high_burden = Float32[]
    # high_burden_sac = Float32[]
    # adult_prev = Float32[]
    #
    # for i in 1 : length(r)
    #     push!(times, r[i].time)
    #     push!(prev, r[i].pop_prev)
    #     push!(sac_prev, r[i].sac_prev)
    #     push!(high_burden, r[i].population_burden[3])
    #     push!(high_burden_sac, r[i].sac_burden[3])
    #     push!(adult_prev, r[i].adult_prev)
    # end
    times = (p->p.time).(r)
    prev = (p->p.pop_prev).(r)
    sac_prev = (p->p.sac_prev).(r)
    high_burden = (p->p.population_burden[3]).(r)
    high_burden_sac = (p->p.sac_burden[3]).(r)
    adult_prev = (p->p.adult_prev).(r)

    plot(times, sac_prev, label  = "SAC prevalence", line=(:black, 0.5, 6, :solid))
    plot!(times, high_burden_sac, label  = "SAC high burden", line=(:purple, 0.5, 6))
    plot!(


        size=(800, 600),

        xticks = (0:100:300),
        yticks = 0:10:100,

        ylabel = "Prevalence",
        xlabel = "Year",


        xrotation = rad2deg(pi/3),

        fillrange = 0,
        fillalpha = 0.25,
        fillcolor = :lightgoldenrod,

        background_color = :ivory
        )
    ylims!((0, 100))

end




function plot_sac_burden_2_runs(r1, r2)
    times1 = Float32[]
    prev1 = Float32[]
    sac_prev1 = Float32[]
    high_burden1 = Float32[]
    high_burden_sac1 = Float32[]
    adult_prev1 = Float32[]

    for i in 1 : length(r1)
        push!(times1, r1[i].time)
        push!(prev1, r1[i].pop_prev)
        push!(sac_prev1, r1[i].sac_prev)
        push!(high_burden1, r1[i].population_burden[3])
        push!(high_burden_sac1, r1[i].sac_burden[3])
        push!(adult_prev1, r1[i].adult_prev)
    end

    times2 = Float32[]
    prev2 = Float32[]
    sac_prev2 = Float32[]
    high_burden2 = Float32[]
    high_burden_sac2 = Float32[]
    adult_prev2 = Float32[]

    for i in 1 : length(r2)
        push!(times2, r2[i].time)
        push!(prev2, r2[i].pop_prev)
        push!(sac_prev2, r2[i].sac_prev)
        push!(high_burden2, r2[i].population_burden[3])
        push!(high_burden_sac2, r2[i].sac_burden[3])
        push!(adult_prev2, r2[i].adult_prev)
    end


    plot(times1, sac_prev1, label  = "SAC prevalence 1", line=(:black, 0.5, 6, :solid))
    plot!(times2, sac_prev2, label  = "SAC prevalence 2", line=(:purple, 0.5, 6))
    plot!(


        size=(800, 600),

        xticks = (0:100:300),
        yticks = 0:10:100,

        ylabel = "Prevalence",
        xlabel = "Year",


        xrotation = rad2deg(pi/3),

        fillrange = 0,
        fillalpha = 0.25,
        fillcolor = :lightgoldenrod,

        background_color = :ivory
        )
    ylims!((0, 100))

end





function get_mean_eggs_age(age_bins, sim_ages, sim_eggs)
    a = Float32[]
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
    a = Float32[]
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
