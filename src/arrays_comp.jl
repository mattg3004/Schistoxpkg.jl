function human_cercariae_maturity(human_cercariae, female_worms, male_worms, time_step)

#=  loop over humans  =#
    for i in 1:size(human_cercariae)[1]

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_cercariae  =#
        if length(human_cercariae[i]) > round(35/time_step, digits = 0)
            females = rand(Binomial(human_cercariae[i][1], 0.5))[1]
            female_worms[i][1] += females
            male_worms[i][1] += (human_cercariae[i][1] - females)
            splice!(human_cercariae[i],1)
        end
    end

#= return arrays  =#
    return human_cercariae, female_worms, male_worms
end


function calculate_worm_pairs(female_worms, male_worms)
    return min.(sum.(female_worms), sum.(male_worms))
end

# function to calculate the number of male and female worms in each human


function calculate_total_worms(female_worms, male_worms)
    return sum.(female_worms), sum.(male_worms)
end


function egg_production(eggs, max_fecundity, r, worm_pairs,
                        total_female_worms, total_male_worms,
                        density_dependent_fecundity, time_step)


# loop over individuals
    for i in 1 : length(eggs)

#= if we have a positive number of worms, then make calculation,
otherwise the number of eggs is trivially 0 =#
        @inbounds if worm_pairs[i] > 0

# calculate the mean number of eggs we would expect
                    @inbounds    mean_eggs =  max_fecundity * worm_pairs[i] *
                        exp(- density_dependent_fecundity *
                        (total_female_worms[i] ))

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


function miracidia_production(eggs, env_miracidia, time_step, age_contact_rate)
#= as we can step forward an arbitrary number of days at a time, we multiply the number of miracidia by the
    length of the forward step, assuming that each of the last given number of days were equivalent to each other
=#
    max_contact_rate = maximum(age_contact_rate)
    xx = age_contact_rate ./ max_contact_rate
    released_eggs = xx .* eggs
    push!(env_miracidia,  sum(released_eggs))
    return env_miracidia
end


function cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate,
    predisposition, age_contact_rate, vac_status, vaccine_effectiveness, human_cercariae_prop)

#= we want the human population to randomly pick up larvae.
therefore we want to shuffle the population.
the following lines make a random permutation of indices for the population=#
    k = size(human_cercariae)[1]
    x = randperm(k)
    uptakes = 0
#= assign larvae which have been in the environment for 40 days to become infective.
then delete those larvae from the environmental larvae =#

    if  length(env_miracidia) > (40 / time_step)
        env_cercariae += (env_miracidia[1] * human_cercariae_prop)
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
                    env_cercariae * time_step / k

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
