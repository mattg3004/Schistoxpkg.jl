

##########################
# function to get age dependent death rate.
# The first entry is for under 1's, and the rest are at 5 year intervals
#  At some point this will be changed to be read from a file

# @param {age_death_rate_per_1000} - death rate per 1000 humans per year.
# Hard coded currently
# @param {time step} - how many days we step forward each simulation time step
# @output{death_rate_per_time_step} - age dependent probability of death each time step

function make_death_rate_array(age_death_rate_per_1000, time_step)

# convert the death rate per 1000 to deaths per day

    death_rate_per_time_step = 1 .- exp.(-1 .* time_step * age_death_rate_per_1000./(1000*365))
# return death rate data
    return death_rate_per_time_step

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

function make_age_contact_rate_array(max_age,  scenario)
    if max_age < 60
        error("max_age must be greater than 60")
    else
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

        return contact_rates_by_age
    end
# return contact rates for age
end









# define a function which will create the initial population.
# This randomly chooses age, and gender.
# predisposition is taken to be gamma distributed.
# There is also a male and female adjustment to predisposition
# adjusting for gender specific behaviour

function create_population(N, max_age, initial_worms, contact_rates_by_age,
    death_rate_per_time_step,worm_stages, female_factor, male_factor,
    initial_miracidia, initial_miracidia_days, predis_aggregation, time_step,
    mda_adherence)

# initialize all the arrays we will keep track of over time

    female_worms = []
    male_worms = []
    human_cercariae = []
    eggs = []
    vac_status = []
    treated = []
    vaccinated = []
    age_contact_rate = []
    death_rate = []
    ages = []
    gender = []
    adherence = []
#=  initialize the Gamma distribution for predisposition selection  =#
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

#=  initialize and fill the environmental variable  =#
    env_miracidia = []

    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia )
    end

    predisposition = rand(gamma_pre,N)

    for i in 1:N

#=  begin pushing entries to the data variables we keep track of  =#
        push!(ages, rand()*max_age)
        push!(gender, rand([0,1]))
        #push!(predisposition, rand(gamma_pre)[1])
        push!(human_cercariae,[])
        push!(eggs,0)
        push!(vac_status, 0)
        push!(treated,0)
        push!(vaccinated, 0)
        if rand() > mda_adherence
            push!(adherence, 0)
        else
            push!(adherence, 1)
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

#=  death rate for the correct age is pushed into array  =#
        push!(death_rate, find_death_rate(age, death_rate_per_time_step))

#=  if the person is chosen to be a male of female, then
    adjust their predisposition based on the
    male or female factor, adjusting for
    behavioural differences according to gender  =#
        x = gender[i]
        predisposition[i] = predisposition[i] * (1-x) *female_factor + predisposition[i] * x *female_factor
    end

#=  return all data that we will use to track the spread of the disease  =#
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status,
            treated, female_worms, male_worms, vaccinated,
            age_contact_rate, death_rate, env_miracidia, adherence
end







#= function to update the contact rate of individuals in the population =#

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
            @inbounds death_rate[i] = death_rate_per_time_step[1]
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

function cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
    predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= we want the human population to randomly pick up larvae.
therefore we want to shuffle the population.
the following lines make a random permutation of indices for the population=#
    k = size(human_cercariae)[1]
    x = randperm(k)
    uptakes = 0
#= assign larvae which have been in the environment for 40 days to become infective.
then delete those larvae from the environmental larvae =#

    if  length(env_miracidia) > (40 / time_step)
        env_cercariae += env_miracidia[1]
        splice!(env_miracidia, 1)
    end

#= loop over the population uptaking larvae=#
    for i in 1:k

# set index to the correct value from the random permutation
        @inbounds  j = x[i]

#= if there are still infective larvae in the environment,
we will uptake from choose from Poisson
distribution. otherwise, just uptake 0. =#
        if env_cercariae > 0
# calculate the rate of the poisson distribution
            @inbounds  pois_rate  = predisposition[j] * contact_rate * age_contact_rate[j] *
                    env_cercariae * time_step

# reduce the rate according to the effectiveness of the vaccine (if any is given)
            @inbounds  pois_rate = pois_rate * (1 - (vac_status[j] > 0) * vaccine_effectiveness);

# choose from the Poisson distribution
            uptake = rand(Poisson(pois_rate))
    #println(uptake)
        else
            uptake = 0
        end

# push the number of uptaken larvae into the human larvae array
        push!(human_cercariae[j], uptake)

# # reduce the infective larvae by the number of larvae uptaken
        env_cercariae -= uptake

        end
#println(uptakes/k)
# return the infective, human and environmental larvae arrays
    return env_cercariae, human_cercariae, env_miracidia

end





# function to kill miracidia in the environment
function miracidia_death(env_miracidia, env_miracidia_death_rate)
    #= as env_miracidia is an array, we need to use . syntax to apply
    the functions to each element in the array =#
    # return rand.(Binomial.(env_miracidia, 1 - env_miracidia_death_rate))
    env_miracidia[end] = trunc(Int, round(env_miracidia[end]/2, digits = 0))
    return env_miracidia
end




# function to kill cercariae in the environment
function cercariae_death(env_cercariae, env_cercariae_death_rate, time_step)
    # updated_cercariae = 0
    # for i in 1:time_step
    #     updated_cercariae += (env_cercariae/time_step) * (1/(2^i))
    # end
    updated_cercariae = trunc(Int,round(env_cercariae / 2, digits= 0))
    return updated_cercariae
end






# Worms die have a specified mean life span, and hence a rate of deaths per day .
# the number of deaths, or maturing from each stage is dependent
# on the number of worm stages and rate of aging through stages and dying.
# This p is multiplied by the time scale of the simulation, so if 2 days
# pass between consecutive time points, twice as many worms age and die


function worm_maturity(female_worms, male_worms, worm_stages,
    average_worm_lifespan, time_step)

# loop over the worms
    for i in 1:size(female_worms)[1]

# probability of aging out of category/ dying
        p = time_step / ( worm_stages * 365 * average_worm_lifespan)

# kill appropriate number of worms in the final stage
        @inbounds n = female_worms[i][worm_stages]
        if(n>0)
            dis = Binomial(n, 1-p)
            @inbounds female_worms[i][worm_stages] = rand(dis, 1)[1]
        end


        @inbounds n = male_worms[i][worm_stages]
        if n > 0
            dis = Binomial(n, 1-p)
            @inbounds male_worms[i][worm_stages] = rand(dis, 1)[1]
        end

#=
     for aging worms, we do this in reverse order, which ensures the
     correct order of aging is respected
=#

        for j in (worm_stages-1):-1:1
#=   choose the number of male and female worms to age from one stage to the next   =#
            @inbounds aging_females = rand(Binomial(female_worms[i][j], p), 1)[1]
            @inbounds aging_males = rand(Binomial(male_worms[i][j], p), 1)[1]

#=   add and subtract the number of worms from the appropriate categories   =#
            @inbounds female_worms[i][j+1] += aging_females
            @inbounds female_worms[i][j] -= aging_females
            @inbounds male_worms[i][j+1] += aging_males
            @inbounds male_worms[i][j] -= aging_males
        end
    end

#=  return the female and male worm arrays  =#
    return female_worms, male_worms

end



# function to calculate the number of male-female pairs of worms in each human
# this is just the minimum of total female worms and the total male worms


function calculate_worm_pairs(female_worms, male_worms)
    return min.(sum.(female_worms), sum.(male_worms))
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
                        total_female_worms, total_male_worms,
                        density_dependent_fecundity, time_step)


# loop over individuals
    for i in 1 : size(worm_pairs)[1]

#= if we have a positive number of worms, then make calculation,
otherwise the number of eggs is trivially 0 =#
        @inbounds if worm_pairs[i] > 0

# calculate the mean number of eggs we would expect
                @inbounds    mean_eggs =  max_fecundity * worm_pairs[i] *
                    exp(- density_dependent_fecundity *
                    (total_female_worms[i] + total_male_worms[i]))

# calculate the number of successes
                @inbounds      NB_r = r * worm_pairs[i]

# calculate the probability of a success
                p = NB_r/(NB_r+mean_eggs)

# choose from NB
                eggs_num = rand(NegativeBinomial(NB_r,p))[1]
            #eggs_num = round(mean_eggs)
            #println("prop = ", eggs_num/mean_eggs)
            else
                eggs_num = 0
            end

# put this selected number of eggs into the eggs array
            @inbounds eggs[i] = eggs_num
        end

# return the eggs array
    return eggs
end


# hatch the eggs in the humans into the environment

function miracidia_production(eggs, env_miracidia, time_step)
#= as we can step forward an arbitrary number of days at a time, we multiply the number of miracidia by the
    length of the forward step, assuming that each of the last given number of days were equivalent to each other
=#
    push!(env_miracidia,  sum(eggs))
    return env_miracidia
end



# function to kill humans at age dependent rate


function death_of_human(ages, gender, predisposition,  human_cercariae, eggs,
                            vac_status, treated, female_worms, male_worms,
                            vaccinated, age_contact_rate, death_rate, time_step,
                            adherence)

#= loop through the population, and based on the death rate,
delete individuals from the population  =#

    for i in size(ages)[1]:-1:1
        r = rand()
        if ages[i] > 100


            if r < time_step/365
                splice!(ages, i)
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
                splice!(death_rate, i)
                splice!(adherence, i)
            end
            else
#= if random number is smaller than the death rate, then delete
    individual from all the arrays  =#
            if r < death_rate[i]
                splice!(ages, i)
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
                splice!(death_rate, i)
                splice!(adherence, i)
            end
      end
    end

#=  return the arrays  =#
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status,
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, adherence
end







# function to add a person to the data if a birth occurs

function birth_of_human(ages, gender, predisposition, human_cercariae, eggs, vac_status,
                        treated, female_worms, male_worms,vaccinated, age_contact_rate,
                        death_rate, female_factor, male_factor, contact_rates_by_age,
                        death_rate_per_time_step, worm_stages, predis_aggregation, adherence,
                        mda_adherence)


#  load the gamma distribution for the predispostion distribution
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)


# fill female and male worms array
    f_worms = fill(0, worm_stages)
    m_worms = fill(0, worm_stages)

# push into arrays
    push!(female_worms, f_worms)
    push!(male_worms, m_worms)
    push!(ages, 0)
    push!(gender, rand([0,1]))
    push!(predisposition, rand(gamma_pre)[1])
    push!(human_cercariae,[])
    push!(eggs,0)
    push!(vac_status, 0)
    push!(treated, 0)
    push!(vaccinated, 0)
    predisp = rand(gamma_pre)[1]

    push!(death_rate, death_rate_per_time_step[1])
    push!(age_contact_rate, contact_rates_by_age[1])
    if rand() > mda_adherence
        push!(adherence, 0)
    else
        push!(adherence, 1)
    end
# adjust predisposition based on gender specific behaviour parameters
    if gender[end] === 0
        predisposition[end] = predisp * female_factor
    else
        predisposition[end] = predisp * male_factor
    end

#  return arrays
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status,
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, adherence
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
    treated, mda_round, gender, adherence)

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
        female_worms = administer_drug(female_worms, y, mda_effectiveness, adherence)
        male_worms = administer_drug(male_worms, y, mda_effectiveness, adherence)
        human_cercariae = administer_drug(human_cercariae, y, mda_effectiveness, adherence)
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
    if mda_round === size(mda_info)[1]
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
    while mda_time < last_mda_time
        push!(mda_info, mda_information(pre_SAC_prop, 0, 4, pre_SAC_gender, mda_effectiveness, mda_time))
        push!(mda_info, mda_information(SAC_prop, 4, 16, SAC_gender, mda_effectiveness, mda_time))
        push!(mda_info, mda_information(adult_prop, 16, 110, adult_gender, mda_effectiveness, mda_time))
        mda_time += regularity
    end
    return mda_info
end



# function to add vaccination to population

function vaccinate(vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
    vaccine_gender, ages, female_worms, male_worms, human_cercariae, eggs,
    treated, vaccine_duration, vac_status, vaccine_round, gender, adherence)

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
        female_worms = administer_drug(female_worms, y, vaccine_effectiveness, adherence)
        male_worms = administer_drug(male_worms, y, vaccine_effectiveness, adherence)
        human_cercariae = administer_drug(human_cercariae, y, vaccine_effectiveness, adherence)
        eggs = administer_drug(eggs, y, 1, adherence)
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


# function to update the variable arrays a given number of times (num_time_steps)

function update_env(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, adherence, mda_adherence,
    record_frequency)

    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = record_frequency
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
            a = get_prevalences(ages, eggs, gamma_k, sim_time)
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
                               total_female_worms, total_male_worms,
                               density_dependent_fecundity, time_step)

#=  mature worms in each human  =#
        female_worms, male_worms = worm_maturity(female_worms, male_worms,
                                                 worm_stages, average_worm_lifespan,
                                                 time_step)

 #=  reduce the vaccination status by the time step  =#
        vac_status = vac_status .- time_step/365

#=  hacth the human eggs into the environment  =#
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step)

#=  update population due to death  =#
       ages , gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate, adherence =
             death_of_human(ages, gender, predisposition,  human_cercariae, eggs,
                                         vac_status, treated, female_worms, male_worms,
                                         vaccinated, age_contact_rate, death_rate, time_step,
                                         adherence)

#=  uptake larvae into humans from the environment  =#
        env_cercariae, human_cercariae, env_miracidia =
                 cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            female_worms, male_worms, human_cercariae, eggs =
            mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender,
                     ages, female_worms, male_worms, human_cercariae, eggs,
                     treated, mda_round, gender, adherence)

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
                    treated, vaccine_duration, vac_status, vaccine_round, gender, adherence)

#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                update_vaccine(vaccine_info, vaccine_round)
        end

#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_death_rate)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_death_rate, time_step)

 #=  choose from binomial distribution for the number of births in the population  =#
        l = rand(Binomial(size(ages)[1], birth_rate))[1]

 #=  loop over the number of births that there are  =#
        if l > 0
            for i in 1:l

 #=  update the population due to births  =#
                 ages , gender, predisposition,  human_cercariae, eggs, vac_status,
                 treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate, adherence =
                     birth_of_human(ages , gender, predisposition,  human_cercariae, eggs, vac_status,
                                        treated, female_worms, male_worms,vaccinated, age_contact_rate,
                                        death_rate, female_factor, male_factor, contact_rates_by_age,
                                        death_rate_per_time_step, worm_stages, predis_aggregation, adherence,
                                        mda_adherence)
            end
        end
    end

#=  return the arrays  =#
    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,
    record
end






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
    pop_prev
    sac_prev
    sac_pop
    final_ages
    time
end


function get_prevalences(ages, eggs, gamma_k, time)

    pop_burden = [0,0,0]
    sac_burden = [0,0,0]
    pop_prev = 0
    sac_prev = 0
    sac_pop = 0
    final_ages = []

    num_humans = size(ages)[1]

    for i in 1:num_humans
        push!(final_ages, ages[i]);
        final_eggs = kato_katz(eggs[i], gamma_k)
        if ages[i] > 4 && ages[i] < 17
            sac_pop = sac_pop + 1;
        end
        if final_eggs > 16
            pop_burden[3] = pop_burden[3] + 1
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if 4 < ages[i] && ages[i] < 17
                sac_burden[3] = sac_burden[3] + 1
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
        elseif final_eggs > 4
            pop_burden[2] = pop_burden[2] + 1
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if ages[i] > 4 && ages[i] < 17
                sac_burden[2] = sac_burden[2] + 1
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev + 1
            end
        elseif final_eggs > 0
            pop_burden[1] = pop_burden[1] + 1
            pop_prev = pop_prev + 1
            if ages[i] > 4 && ages[i] < 17
                sac_burden[1] = sac_burden[1] + 1
                sac_prev = sac_prev +  1
            end
        end
    end

    output = out( round.(100 .*pop_burden./num_humans, digits = 2),
        round.(100 .*sac_burden./sac_pop, digits = 2),
        round.(100 .*pop_prev ./ num_humans, digits = 2),
        round(100 .*sac_prev / sac_pop, digits = 2),
        sac_pop, round(mean(final_ages), digits = 2),
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
    death_rate = get!(d, "death_rate",1)
    env_miracidia = get!(d, "env_miracidia",1)
    env_cercariae = get!(d, "env_cercariae",1)
    adherence = get!(d, "adherence",1)

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
        death_rate = death_rate[x]
        adherence = adherence[x]
    end
     return ages, gender, predisposition, human_cercariae, eggs, vac_status, treated,
        female_worms, male_worms, vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence
end



function save_population_to_file(filename, ages, gender, predisposition, human_cercariae, eggs, vac_status, treated,
        female_worms, male_worms, vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae)
    save(filename, "ages", ages ,  "gender", gender,"predisposition",   predisposition,
        "human_cercariae", human_cercariae, "eggs", eggs,
        "vac_status", vac_status,"treated", treated, "female_worms",  female_worms, "male_worms", male_worms,
        "vaccinated", vaccinated,  "age_contact_rate", age_contact_rate, "death_rate", death_rate,
        "env_miracidia",env_miracidia, "env_cercariae", env_cercariae, "adherence", adherence)
end



function run_simulation_from_loaded_population(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info,
    record_frequency)


    ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms, vaccinated,
    age_contact_rate, death_rate, env_miracidia, env_cercariae, record =
    update_env(num_time_steps, ages, human_cercariae, female_worms, male_worms,
        time_step, average_worm_lifespan,
        eggs, max_fecundity, r, worm_stages,
        vac_status, gender, predis_aggregation,
        predisposition, treated, vaccine_effectiveness,
        density_dependent_fecundity,
        vaccinated, age_contact_rate, death_rate, env_miracidia,
        env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
        female_factor, male_factor, contact_rates_by_age,
        death_rate_per_time_step, birth_rate, mda_info, vaccine_info, adherence, mda_adherence,
        record_frequency)

    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms, vaccinated,
    age_contact_rate, death_rate, env_miracidia, env_cercariae, record, adherence

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
    env_cercariae_death_rate, env_miracidia_death_rate, mda_coverage, mda_round, vaccine_effectiveness,
    mda_info, vaccine_info, record_frequency, mda_adherence, scenario)


    age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                               0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                               21.83, 29.98, 36.98]


    contact_rates_by_age = make_age_contact_rate_array(max_age, scenario)
    death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)


    ages , gender, predisposition,  human_cercariae, eggs, vac_status,
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_miracidia, adherence =
        create_population(N, max_age, initial_worms, contact_rates_by_age,
        death_rate_per_time_step, worm_stages, female_factor, male_factor,
        initial_miracidia, initial_miracidia_days, predis_aggregation, time_step,
        mda_adherence)



    ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms, vaccinated,
    age_contact_rate, death_rate, env_miracidia, env_cercariae, record =
        update_env(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, adherence,
    record_frequency)

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
    ages = []
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
function create_population_specified_ages(N, initial_worms, contact_rates_by_age,
        worm_stages, female_factor, male_factor,initial_miracidia,
        initial_miracidia_days, predis_aggregation, time_step,
        spec_ages, ages_per_index, death_rate_per_time_step,
        mda_adherence)

# initialize all the arrays we will keep track of over time

    female_worms = []
    male_worms = []
    human_cercariae = []
    eggs = []
    vac_status = []
    treated = []
    vaccinated = []
    age_contact_rate = []
    death_rate = []
    gender = []
    adherence = []
#=  initialize the Gamma distribution for predisposition selection  =#
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

#=  initialize and fill the environmental variable  =#
    env_miracidia = []

    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia )
    end

    predisposition = rand(gamma_pre,N)
    ages = specified_age_distribution(N, spec_ages, ages_per_index)

    for i in 1:N

#=  begin pushing entries to the data variables we keep track of  =#
        push!(gender, rand([0,1]))
        push!(human_cercariae,[])
        push!(eggs,0)
        push!(vac_status, 0)
        push!(treated,0)
        push!(vaccinated, 0)
        if rand() > mda_adherence
            push!(adherence, 0)
        else
            push!(adherence, 1)
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

#=  death rate for the correct age is pushed into array  =#
        push!(death_rate, find_death_rate(age, death_rate_per_time_step))

#=  if the person is chosen to be a male of female, then
    adjust their predisposition based on the
    male or female factor, adjusting for
    behavioural differences according to gender  =#
        x = gender[i]
        predisposition[i] = predisposition[i] * (1-x) *female_factor + predisposition[i] * x *female_factor
    end

#=  return all data that we will use to track the spread of the disease  =#
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status,
            treated, female_worms, male_worms, death_rate, age_contact_rate,
            vaccinated, env_miracidia, adherence
end



#= this function will run the population for a specified number of time steps
with no births, deaths or aging in the population, hence running it to equilibrium
for the specified size of population =#

function update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,vaccinated, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate)


    sim_time = 0
    record_time = record_frequency
    record = []


#=  loop for number of sims  =#

    for j in 1:num_time_steps


        if sim_time >= record_time
            a = get_prevalences(ages, eggs, gamma_k, sim_time)
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



#=  hacth the human eggs into the environment  =#
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step)

#=  uptake larvae into humans from the environment  =#
        env_cercariae, human_cercariae, env_miracidia =
                 cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)


#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_death_rate)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_death_rate, time_step)


    end

#=  return the arrays  =#
    return ages, gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, env_miracidia, env_cercariae, record
end





# function to update the variable arrays a given number of times (num_time_steps) where no aging, births or
# deaths occur. This differs from update_env_to_equlibrium as mda and vaccination can take place


function update_env_no_births_deaths(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, adherence, mda_adherence,
    record_frequency)

    update_contact_death_rates = 1/5
    sim_time = 0
    record_time = record_frequency
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
            a = get_prevalences(ages, eggs, gamma_k, sim_time)
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
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step)


#=  uptake larvae into humans from the environment  =#
        env_cercariae, human_cercariae, env_miracidia =
                 cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            female_worms, male_worms, human_cercariae, eggs =
            mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender,
                     ages, female_worms, male_worms, human_cercariae, eggs,
                     treated, mda_round, gender, adherence)

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
                    treated, vaccine_duration, vac_status, vaccine_round, gender, adherence)

#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                update_vaccine(vaccine_info, vaccine_round)
        end

#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_death_rate)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_death_rate, time_step)


    end

#=  return the arrays  =#
    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,
    record
end




# allow mda and vaccine when updating the model, but whenever there is a death, we instantly add another individual
function update_env_keep_population_same(num_time_steps, ages, human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, adherence, mda_adherence,
    record_frequency)

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
            a = get_prevalences(ages, eggs, gamma_k, sim_time)
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
        env_miracidia = miracidia_production(eggs, env_miracidia, time_step)

    #=  update population due to death  =#
       ages , gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate, adherence =
             death_of_human(ages, gender, predisposition,  human_cercariae, eggs,
                                         vac_status, treated, female_worms, male_worms,
                                         vaccinated, age_contact_rate, death_rate, time_step,
                                         adherence)

        if length(ages) < start_pop
            num_births = start_pop - length(ages)
            for i in 1:num_births
                ages , gender, predisposition,  human_cercariae, eggs, vac_status,
                 treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate, adherence =
                     birth_of_human(ages , gender, predisposition,  human_cercariae, eggs, vac_status,
                                        treated, female_worms, male_worms,vaccinated, age_contact_rate,
                                        death_rate, female_factor, male_factor, contact_rates_by_age,
                                        death_rate_per_time_step, worm_stages, predis_aggregation, adherence,
                                        mda_adherence)
            end
        end
#=  uptake larvae into humans from the environment  =#
        env_cercariae, human_cercariae, env_miracidia =
                 cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= check if we are at a point in time in which an mda is scheduled to take place =#
        if sim_time >= next_mda_time

#= perform mda =#
            female_worms, male_worms, human_cercariae, eggs =
            mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender,
                     ages, female_worms, male_worms, human_cercariae, eggs,
                     treated, mda_round, gender, adherence)

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
                    treated, vaccine_duration, vac_status, vaccine_round, gender, adherence)

#= update information for the next round of vaccination =#
            vaccine_round += 1
            vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
                                update_vaccine(vaccine_info, vaccine_round)
        end

#=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_death_rate)

#=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_death_rate, time_step)


    end

#=  return the arrays  =#
    return ages , gender, predisposition,  human_cercariae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,
    record
end




# when we run multiple simulations, we store them in an array. This function will store the prevalence and sac prevalence
function collect_prevs(times, prev, sac_prev, record, run)
        if run == 1
            for i in 1 : length(record)
                push!(times, record[i].time)
                push!(prev, [record[i].pop_prev])
                push!(sac_prev, [record[i].sac_prev])
            end
        else
            for i in 1 : length(record)
                push!(prev[i], record[i].pop_prev)
                push!(sac_prev[i], record[i].sac_prev)
            end
        end
    return times, prev, sac_prev
end

# repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death
function run_repeated_sims_no_population_change(num_repeats, num_time_steps,
    time_step, average_worm_lifespan,
    max_fecundity, r, worm_stages, predis_aggregation, vaccine_effectiveness,
    density_dependent_fecundity, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age,
    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, mda_adherence,
    record_frequency, times, prev, sac_prev, filename)




    for run in 1:num_repeats


        ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, death_rate_equ, env_miracidia_equ , env_cercariae_equ, adherence_equ =
        load_population_from_file(filename, N, true)

        ages , gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,
        record =
         update_env_keep_population_same(num_time_steps, copy(ages_equ), copy(human_cercariae_equ), copy(female_worms_equ), copy(male_worms_equ),
                    time_step, average_worm_lifespan,
                    copy(eggs_equ), max_fecundity, r, worm_stages,
                    copy(vac_status_equ), copy(gender_equ), predis_aggregation,
                    copy(predisposition_equ), copy(treated_equ), vaccine_effectiveness,
                    density_dependent_fecundity,
                    copy(vaccinated_equ) , copy(age_contact_rate_equ), copy(death_rate_equ), copy(env_miracidia_equ) ,
                    copy(env_cercariae_equ) , contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
                    female_factor, male_factor, contact_rates_by_age,
                    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, copy(adherence_equ), mda_adherence,
                    record_frequency);

        times, prev, sac_prev = collect_prevs(times, prev, sac_prev, record, run)
        if mod(run,10)==0
            println("Done ", run)
        end

    end
    return times, prev, sac_prev
end




# repeat simulations where we allow mdas and vaccination, and allow random births and deaths
function run_repeated_sims_random_births_deaths(num_repeats, num_time_steps,
                time_step, average_worm_lifespan,
                max_fecundity, r, worm_stages,
                predis_aggregation,
                density_dependent_fecundity,
                age_contact_rate_equ, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
                female_factor, male_factor, contact_rates_by_age,
                death_rate_per_time_step, birth_rate, mda_info, vaccine_info,  mda_adherence,
                record_frequency, times, prev, sac_prev, filename)




    for run in 1:num_repeats

        ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, death_rate_equ, env_miracidia_equ , env_cercariae_equ, adherence_equ =
        load_population_from_file(filename, N, true)

        ages , gender, predisposition,  human_cercariae, eggs,
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae, adherence,
        record =
         update_env(num_time_steps, copy(ages_equ), copy(human_cercariae_equ), copy(female_worms_equ), copy(male_worms_equ),
                    time_step, average_worm_lifespan,
                    copy(eggs_equ), max_fecundity, r, worm_stages,
                    copy(vac_status_equ), copy(gender_equ), predis_aggregation,
                    copy(predisposition_equ), copy(treated_equ), vaccine_effectiveness,
                    density_dependent_fecundity,
                    copy(vaccinated_equ) , copy(age_contact_rate_equ), copy(death_rate_equ), copy(env_miracidia_equ) ,
                    copy(env_cercariae_equ) , contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
                    female_factor, male_factor, contact_rates_by_age,
                    death_rate_per_time_step, birth_rate, mda_info, vaccine_info, copy(adherence_equ), mda_adherence,
                    record_frequency);


        times, prev, sac_prev = collect_prevs(times, prev, sac_prev, record, run)
        if mod(run,10)==0
            println("Done ", run)
        end

    end
    return times, prev, sac_prev
end
