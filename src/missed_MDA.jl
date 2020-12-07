using Schistoxpkg
using Random
using JLD
using Distributions
using PyPlot


function plot_sac_burden_and_sac_high_burden(rec, years)
    times = []
    prev = []
    sac_prev = []
    high_burden = []
    high_burden_sac = []
    adult_prev = []

    for i in 1 : length(rec)
        push!(times, rec[i].time)
        push!(prev, rec[i].pop_prev)
        push!(sac_prev, rec[i].sac_prev)
        push!(high_burden, rec[i].population_burden[3])
        push!(high_burden_sac, rec[i].sac_burden[3])
        push!(adult_prev, rec[i].adult_prev)
    end

    plt.plot(times, sac_prev)
    # plot!(times, high_burden_sac, label  = "SAC high burden", line=(:purple, 0.5, 6))
    # plot!(
    #
    #
    #     size=(800, 600),
    #
    #     xticks = (0:100:years),
    #     yticks = 0:10:100,
    #
    #     ylabel = "Prevalence",
    #     xlabel = "Year",
    #
    #
    #     xrotation = rad2deg(pi/3),
    #
    #     fillrange = 0,
    #     fillalpha = 0.25,
    #     fillcolor = :lightgoldenrod,
    #
    #     background_color = :ivory
    #     )
    # ylims!((0, 100))
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

filename = "high_prev_high_adult_missed_mdas.jld"
scenario = "moderate adult" # must be one of "high adult", "low adult" and "moderate adult"


fig1 = "high_adult_1_MDA.png"
fig2 = "high_adult_SAC_1_MDA.png"
fig3 = "high_adult_high_burden_1_MDA.png"
fig4 = "high_adult_high_burden_sac_1_MDA.png"
fig5 = "high_adult_5_MDA.png"
fig6 = "high_adult_SAC_5_MDA.png"
fig7 = "high_adult_high_burden_5_MDA.png"
fig8 = "high_adult_high_burden_sac_5_MDA.png"

include("parameters.jl")


num_repeats = 10

predis_aggregation = 0.24
predis_weight = 1
 # drug_efficacy = 1

r = 0.035 # aggregation parameter for negative binomial for egg production
# drug_efficacy = 1

number_years_equ = 100

num_time_steps_equ = trunc(Int, 365*number_years_equ / time_step)
contact_rate = 0.2
human_cercariae_prop = 1
env_cercariae_survival_prop = 1/2
env_miracidia_survival_prop = 1/2
N = 1000
max_fecundity = 0.24
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [], []);

contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, [1,2,3], [3,5,6]);
N_communities = 1
community_probs = 1
initial_miracidia = 50000*N/1000
init_env_cercariae = 50000*N/1000
community_contact_rate = 1
worm_stages=2
miracidia_maturity_time = 24
## run to equilibrium
ages , death_ages, gender, predisposition, community,
   human_cercariae, eggs, vac_status,
   treated, female_worms, male_worms, age_contact_rate,
   vaccinated, env_miracidia, adherence, access = create_population_specified_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                worm_stages, female_factor, male_factor,initial_miracidia,
                initial_miracidia_days, predis_aggregation, predis_weight,
                time_step,
                spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                mda_adherence, mda_access)

                    death_ages = []
                    for i in 1:N
                         push!(death_ages, get_death_age(death_prob_by_age, ages_for_deaths))
                    end
                    mean(death_ages)
                    ages, death_ages = generate_ages_and_deaths(20000, ages, death_ages, death_prob_by_age, ages_for_deaths)
age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)

mda_info = []
vaccine_info = []
ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ, human_cercariae_equ, eggs_equ,
vac_status_equ, treated_equ, female_worms_equ, male_worms_equ,
vaccinated_equ, age_contact_rate_equ,
env_miracidia_equ, env_cercariae_equ, adherence_equ,access_equ,
record = update_env_keep_population_same(7600, copy(ages), copy(death_ages), copy(community), community_contact_rate, community_probs,
    copy(human_cercariae), copy(female_worms), copy(male_worms),
    time_step, average_worm_lifespan,
    copy(eggs), max_fecundity, r, worm_stages,
    copy(vac_status), copy(gender), predis_aggregation,predis_weight,
    copy(predisposition), copy(treated), vaccine_effectiveness,
    density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
    copy(vaccinated), copy(age_contact_rate), copy(env_miracidia),
    copy(env_cercariae), contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
    record_frequency, human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold,
    kato_katz_par, use_kato_katz)

#  @time ages_equ, gender_equ, predisposition_equ,  human_cercariae_equ, eggs_equ,
# vac_status_equ, treated_equ, female_worms_equ, male_worms_equ,
# vaccinated_equ, env_miracidia_equ, env_cercariae_equ, record_high =
#     update_env_to_equilibrium(7600, copy(ages), copy(human_cercariae), copy(female_worms), copy(male_worms),
#     copy(community), community_contact_rate,
#         time_step, average_worm_lifespan,
#         copy(eggs), max_fecundity, r, worm_stages,
#         copy(vac_status), copy(gender), predis_aggregation,
#         copy(predisposition), copy(treated), vaccine_effectiveness,
#         density_dependent_fecundity, copy(vaccinated), copy(env_miracidia),
#         copy(env_cercariae), contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
#         female_factor, male_factor, contact_rates_by_age, record_frequency, copy(age_contact_rate),human_cercariae_prop,
#         miracidia_maturity_time)

times = (p ->p.time).(record)
sac_burden = (p ->p.sac_burden[1]).(record)
high_sac_burden = (p ->p.sac_burden[3]).(record)
clf()
plt.plot(times, sac_burden)
plt.plot(times, high_sac_burden)
gcf()
plot_sac_burden_and_sac_high_burden(record_high, number_years_equ)

        ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ, human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
         env_cercariae_equ, adherence_equ, access_equ =
        load_population_from_file(filename, N, true)


save(filename, "ages", ages_equ ,  "gender", gender_equ,"predisposition",   predisposition_equ,
    "human_cercariae", human_cercariae_equ, "eggs", eggs_equ, "community", community,
    "vac_status", vac_status_equ,"treated", treated_equ, "female_worms",  female_worms_equ, "male_worms", male_worms_equ,
    "vaccinated", vaccinated_equ,  "age_contact_rate", age_contact_rate, "death_ages", death_ages,
    "env_miracidia",env_miracidia_equ, "env_cercariae", env_cercariae_equ, "adherence", adherence, "access", access)


###### run baseline case where mda isn't missed ########

number_years = 20
num_time_steps = trunc(Int, 365*number_years / time_step)

mda_info = create_mda(0, .75, 0, 1,
    number_years, 1, [0,1], [0,1], [0,1], drug_efficacy)

vaccine_info = []

num_repeats = 10


times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden = run_repeated_sims_no_births_deaths(num_repeats, num_time_steps,
    time_step, average_worm_lifespan, community_contact_rate, community_probs,
    max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
    density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
    record_frequency, filename, human_cercariae_prop, miracidia_maturity_time,  heavy_burden_threshold,
    kato_katz_par, use_kato_katz)

mean_prev_baseline = mean.(prev)
mean_sac_prev_baseline = mean.(sac_prev)
mean_high_burden_baseline = mean.(high_burden)
mean_high_burden_sac_baseline = mean.(high_burden_sac)
mean_adult_prev_baseline = mean.(adult_prev)

clf()
plt.plot(times, mean_sac_prev_baseline )
plt.plot(times, mean_high_burden_sac_baseline )
gcf()
Plots.plot(times, [mean_sac_prev_baseline, mean_high_burden_sac_baseline],
            label = [ "SAC" "Heavy Burden SAC"], dpi = 300)
Plots.plot!([1,1], seriestype="hline", label = "Heavy burden goal")

###### run where mda is missed in the second year ########

ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
 eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
 vaccinated_equ, age_contact_rate_equ, death_rate_equ,
 env_miracidia_equ , env_cercariae_equ, adherence_equ, access_equ =
load_population_from_file(filename, N, true)

# when will the mda be paused and restarted
mda_end = 1
mda_restart = 3

mda_info = create_mda(0, .75, 0, 1,
    mda_end, 1, [0,1], [0,1], [0,1], drug_efficacy)


append!(mda_info, create_mda(0, 0.75, 0, mda_restart,
  number_years, 1, [0,1], [0,1], [0,1], drug_efficacy))



  times, prev_1, sac_prev_1, high_burden_1, high_burden_sac_1, adult_prev_1 =
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
      time_step, average_worm_lifespan,community_contact_rate, community_probs,
      max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
      density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
      female_factor, male_factor, contact_rates_by_age,
      death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
      record_frequency, filename, human_cercariae_prop)


mean_prev_1 = mean.(prev_1)
mean_sac_prev_1 = mean.(sac_prev_1)
mean_high_burden_1 = mean.(high_burden_1)
mean_high_burden_sac_1 = mean.(high_burden_sac_1)
mean_adult_prev_1 = mean.(adult_prev_1)

#
# plot(times, [mean_prev_baseline, mean_prev_1], label = ["Baseline" "Interrupted"],dpi = 300)
# plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
# savefig(fig1)
#
plot(times, [mean_sac_prev_baseline, mean_sac_prev_1], label = ["Baseline SAC" "Interrupted SAC"],dpi = 300)
plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
# savefig(fig2)
#
# plot(times, [mean_high_burden_baseline, mean_high_burden_1], label = ["Baseline high burden" "Interrupted high burden"],dpi = 300)
# plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
# savefig(fig3)

plot(times, [mean_high_burden_sac_baseline, mean_high_burden_sac_1],
label = ["Baseline high burden SAC" "Interrupted high burden SAC"], dpi = 300)
plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
plot!([1,1], seriestype="hline", label = "1%")
savefig(fig4)



###### run where mda is missed after 5 years ########

ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
 eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
 vaccinated_equ, age_contact_rate_equ, death_rate_equ, env_miracidia_equ ,
 env_cercariae_equ, adherence_equ, access_equ =
load_population_from_file(filename, N, true)
# when will the mda be paused and restarted
mda_end = 5
mda_restart = 7

mda_info = create_mda(0, .75, 0, 1,
    mda_end, 1, [0,1], [0,1], [0,1], drug_efficacy)


append!(mda_info, create_mda(0, 0.75, 0, mda_restart,
  number_years, 1, [0,1], [0,1], [0,1], drug_efficacy))



  times, prev_5, sac_prev_5, high_burden_5, high_burden_sac_5, adult_prev_5 = run_repeated_sims_no_population_change(num_repeats, num_time_steps,
      time_step, average_worm_lifespan,community_contact_rate, community_probs,
      max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
      density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
      female_factor, male_factor, contact_rates_by_age,
      death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
      record_frequency, filename, human_cercariae_prop, miracidia_maturity_time)


mean_prev_5 = mean.(prev_5)
mean_sac_prev_5 = mean.(sac_prev_5)
mean_high_burden_5 = mean.(high_burden_5)
mean_high_burden_sac_5 = mean.(high_burden_sac_5)
mean_adult_prev_5 = mean.(adult_prev_5)
#
# plot(times, [mean_prev_baseline, mean_prev_5], label = ["Baseline" "Interrupted" ],dpi = 300)
# plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
# ylims!((0, maximum(mean_prev_baseline)))
# savefig(fig5)
#
 plot(times, [mean_sac_prev_baseline, mean_sac_prev_5], label = ["Baseline SAC" "Interrupted SAC"],dpi = 300)
 plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
 ylims!((0, maximum(mean_sac_prev_baseline)))
# savefig(fig6)
#
# plot(times, [mean_high_burden_baseline, mean_high_burden_5], label = ["Baseline high burden" "Interrupted high burden"],dpi = 300)
# plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
# ylims!((0, maximum(mean_high_burden_baseline)))
# savefig(fig7)

plot(times, [mean_high_burden_sac_baseline, mean_high_burden_sac_5],
label = ["Baseline high burden SAC" "Interrupted high burden SAC"], dpi = 300)
plot!([mda_end+1, mda_restart], seriestype="vline", label = "Missed MDA")
plot!([1,1], seriestype="hline", label = "1%")
savefig(fig8)
