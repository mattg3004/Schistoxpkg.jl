using XLSX
xf = XLSX.readxlsx("/Users/matthewgraham/Dropbox/schistox/papers/gurarie_data_additional_file5.xlsx")
XLSX.sheetnames(xf)

sh = xf["Age_and_Egg_count_Milalani_2009"]

rec_eggs = sh["C2:C778"]
rec_ages = sh["B2:B778"]
length(rec_eggs)
fit_to = []
for i in 1:length(rec_eggs)
    push!(fit_to, rec_eggs[i])
end
histogram(fit_to, bins = 200)


maximum(fit_to)
fit_to_round = round.(fit_to, digits = 0)


rec_eggs_dist = []
for i in 0: maximum(fit_to_round)
    x = findall(fit_to_round .== i)
    push!(rec_eggs_dist, length(x)/length(fit_to_round))
end


sim_eggs_dist = []
for i in 0: maximum(fit_to_round)
    x = findall(sim_eggs .== i)
    push!(sim_eggs_dist, length(x)/length(sim_eggs))
end

plot(0: maximum(fit_to_round),[cumsum(sim_eggs_dist) cumsum(rec_eggs_dist)])


#############################################################################
age = 20

x = findall( age .<= sim_ages_low .<= (age+1 - 0.0001))
worm_pairs = calculate_worm_pairs(female_worms_low[x], male_worms_low[x])
mean(worm_pairs)

total_female_worms, total_male_worms =
 calculate_total_worms(female_worms[x], male_worms[x])


new_eggs = egg_production(eggs[x], max_fecundity, r, worm_pairs ,
                         total_female_worms, total_male_worms,
                         density_dependent_fecundity, time_step)
sum(new_eggs .> 0)/length(new_eggs)
sum(new_eggs .> 16)/length(new_eggs)



#############################################################################
sim_ages_low = record_low[end].final_ages
sim_eggs_low = record_low[end].recorded_eggs
#
# sim_ages_low_slow = record_low_slow[end].final_ages
# sim_eggs_low_slow = record_low_slow[end].recorded_eggs

sim_ages_high = record_high[end].final_ages
sim_eggs_high = record_high[end].recorded_eggs

###############################################################################


age_bins = [0, 5, 10, 15,20,25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
age_bins = 0:3:100





c_low = get_mean_worms_age(age_bins, sim_ages_low, male_worms_low, female_worms_low)
# c_low_slow = get_mean_worms_age(age_bins, sim_ages_low_slow, male_worms_low_slow, female_worms_low_slow)
c_high = get_mean_worms_age(age_bins, sim_ages_high, male_worms_equ, female_worms_equ)

Plots.plot(age_bins, [c_high,c_low])
Plots.plot(age_bins, c_high)

a = get_mean_eggs_age(age_bins, sim_ages_low_slow, sim_eggs_low_slow)
plot( a)


a = get_mean_eggs_age(age_bins, sim_ages_high, sim_eggs_high)
plot(age_bins, a)

b = get_mean_eggs_age(age_bins, rec_ages, rec_eggs)
plot( [a, b])
