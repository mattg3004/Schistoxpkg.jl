
contact_rates = make_age_contact_rate_array(100, "high adult", [4,9,15,100], [0.01, 1.2, 0.5, 0.002])
input_contact_rates = [1,2,3,4]
input_ages = [5,16,25,70]

contact_rates_by_age = [fill(input_contact_rates[end], max_age+1)]
contact_rates_by_age = contact_rates_by_age[1]
println("asgfsadf ", contact_rates_by_age)
println("asgfsadf ", contact_rates_by_age[11])
for i in 1 : length(input_contact_rates)
    if i == 1
        for j in 1:(input_ages[i] + 1)
            print("j = ", j)
            println(input_contact_rates[i])
            println("contact rate = ", contact_rates_by_age[j])
            contact_rates_by_age[j] = input_contact_rates[i]
        end
    else
        for j in (input_ages[i-1]+2):(input_ages[i] + 1)
            contact_rates_by_age[j] = input_contact_rates[i]
        end
    end
end
