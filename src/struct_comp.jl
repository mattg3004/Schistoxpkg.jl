

function create_population_specified_ages(N, initial_worms, contact_rates_by_age,
        worm_stages, female_factor, male_factor,initial_miracidia,
        initial_miracidia_days, predis_aggregation, time_step,
        spec_ages, ages_per_index,
        mda_adherence, mda_access)


#=  initialize the Gamma distribution for predisposition selection  =#
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

#=  initialize and fill the environmental variable  =#
    env_miracidia = []

    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia )
    end

    predisposition = rand(gamma_pre,N)
    ages = specified_age_distribution(N, spec_ages, ages_per_index)
    humans = []
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


function calculate_worm_pairs(d)
    return min(sum(d.female_worms), sum(d.male_worms))
end

# function to calculate the number of male and female worms in each human


function calculate_total_worms(female_worms, male_worms)
    return sum(female_worms) + sum(male_worms)
end



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
#         new_miracidia = new_miracidia + d[i].eggs * d[i].age_contact_rate / max_contact_rate
                new_miracidia = new_miracidia + d[i].eggs
    end
    push!(env_miracidia, new_miracidia)
    return env_miracidia
end


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
