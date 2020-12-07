
function human_larvae_maturity2(humans, pars)

    a = "enact_maturity_function_"
#=  loop over humans  =#
    @inbounds for h in humans

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_cercariae  =#
        fn = a * (length(h.larvae) > round(pars.human_larvae_maturity_time/pars.time_step, digits = 0))
        h = getfield(Schistoxpkg, Symbol(fn))(h)
    end

#= return arrays  =#
    return humans
end

humans, miracidia, cercariae = create_population_specified_ages(pars)

h = humans[1]
push!(h.larvae, 100)
push!(h.larvae, 100)
push!(h.larvae, 100)
push!(h.larvae, 100)
push!(h.larvae, 100)

a = "enact_maturity_function_"
fn = a * string(length(h.larvae) > round(pars.human_larvae_maturity_time/pars.time_step))
h = getfield(Schistoxpkg, Symbol(fn))(h)
fn = a * string(length(h.larvae) > round(pars.human_larvae_maturity_time/pars.time_step))
h = getfield(Schistoxpkg, Symbol(fn))(h)


enact_maturity_function_false(10)



x = 1E-10
