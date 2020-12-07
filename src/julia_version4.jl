

module Schisto_Sim


  #using PyPlot           # commented out because PyPlot not properly installed on Magnus's CPU
  using Distributions
  using ProgressMeter
  using DataFrames
  using JLD
  using CSV
  #using BensEconomicFuns  #can't remember what this is


  # Export functions/structs to use in REPL (e.g. in Jupyyer notebook)

  #export run_MDA, run_vaccine, run_to_equilibrium, run_strategies1, run_strategies2, reintroduction
  export load_age_death_rates, test_run, get_distributions, create_population
  export load_contact_data, load_environment, create_programs, set_beta
  export run_vaccination, vaccine_setup, MDA_test, run_vacmda
  export Human, Environment, MyParameters, Vaccination_Program, Statistics, Epi_distributions, Output
  export generate_age_dist, vacmda_setup, create_VACMDA, find_stable_prevalence
  export reintroduction
  export efficacy_programs
  export create_extra_programs
  #export make_high_prevalence_plot

  mutable struct Human
    age::Int64                         # age in days
    sex::Int64                         # sex, not currently used (should probably be binary)
    predisposition::Float64            # risk factor
    male_worms::Array{Int64,1}         # male worms array, each element being a certain life-stage
    female_worms::Array{Int64,1}       # female worms array
    larvae::Array{Int64,1}             # larvae taken in from environment, cant yet reproduce
    eggs::Int64                        # eggs produced on current day
    vac_status::Int64                  # individual vaccine status. 0 = unvaccinted, >1 = vaccinated (this counts down)
    snc_prob::Float64                  # probability from Louise's SNC paper

    # Constructors
    Human() = r(1,0,1,[0,0,0,0],[0,0,0,0],[],0,0)
    Human(age,sex,predisposition,male_worms,female_worms,larvae,eggs,vac_status,snc_prob) =
      new(age,sex,predisposition,male_worms,female_worms,larvae,eggs,vac_status,snc_prob)
  end

  mutable struct Environment
      larvae::Array{Int64,1}
      infective_larvae::Int64         # FOR COMPARISON WITH ICL MODEL
      humans::Array{Human,1}          # container for population
      time::Int64                     # in days

      # Constructors
      Environment() = new([],0,[],0)
      Environment(larvae,infective_larvae,humans,time) = new(larvae,infective_larvae,humans,time)
  end

  mutable struct MyParameters
      worm_stages::Int64                        # THIS IS NEW, FOR COMPARISON WITH ICL MODEL (changes gamma distributed lifespans)
      death_human::Array{Float64,2}             # human death rate per day by age
      death_worm::Float64                       # worm death rate in humans per day
      contact_rate::Float64                     # contact between humans and larvae per day per larvae per human
      max_fecundity::Float64                    # maximum number of eggs produced per mating female worm
      gamma_r::Distributions.Gamma{Float64}     # distribution of lifelong predisposition factor
      gamma_k::Distributions.Gamma{Float64}     # distribution for kato_katz test
      age_contact::Array{Float64,2}             # relative contact rates by age
      snc_correlation::Float64                  # parameters from Louise's SNC paper
      snc_coverage::Float64

      # constructors with default values
      MyParameters() = new(1,Array{Float64, 2}(0,0), 0.0, 0.00001, 0.05,
        Distributions.Gamma(0.04,1.0/0.04), Distributions.Gamma(0.87,1.0/0.87),
        Array{Float64, 2}(0,0) ,0.05 ,0.75)  #0.55~6.1 0.6~2.1

      MyParameters(worm_stages,death_human,death_worm,contact_rate,max_fecundity,
        gamma_r,gamma_k,age_contact,snc_correlation,snc_coverage) =
        new(death_human,death_worm,contact_rate,max_fecundity,
          gamma_r,gamma_k,age_contact,snc_correlation,snc_coverage)
  end

  # container for histograms
  mutable struct Epi_distributions
      age_histogram::DataFrames.DataFrame
      burden_histogram::DataFrames.DataFrame
      egg_histogram::DataFrames.DataFrame
      prevalence_histogram::DataFrames.DataFrame
      n_samples::Int64
  end



  # container for pop-level statistics
  mutable struct Statistics
      mean_burden::Array{Float64,2}
      prevalence::Array{Float64,2}
      aggregation::Array{Float64,2}
      worms_in_pool::Array{Float64,2}
      stratified_prevalence::Array{Float64,3}
      adult_prevalence::Array{Float64,2}
      adult_stratified_prevalence::Array{Float64,3}


  end

  # container for vaccination program parameters
  mutable struct Vaccination_Program
    #coverage::Float64
    duration::Float64
    target_ages::Array{Int64,1}
    coverages::Array{Float64,1}  # target age needs a coverage associated with it
    years::Int64
    contact_eff::Float64
    fecundity_eff::Float64
    death_eff::Float64
    theraputic_eff::Float64
    intensity::String
  end

  mutable struct MDA_Program
    coverage::Array{Float64,1}
    #target_ages::Array{Int64,1}
    years::Int64
    efficacy::Float64
    R0::Float64
    frequency::Int64

    MDA_Program() = new(0.75*ones(100),16,0.86,1.0,1)
    MDA_Program(coverage,years,efficacy,R0,frequency) =
     new(coverage,years,efficacy,R0,frequency)
  end

  mutable struct VACMDA_Program
    duration::Float64
    target_ages::Array{Int64,1}
    coverages::Array{Float64,1}  # target age needs a coverage associated with it
    years::Int64
    contact_eff::Float64
    fecundity_eff::Float64
    death_eff::Float64
    theraputic_eff::Float64
    intensity::String
    MDA_cov::Float64

  end

  mutable struct Output
    stats::Statistics
    pre_eggs::Array{Float64,1}
    post_eggs::Array{Float64,1}
  end


#include("economic_funs.jl")


#=================================================================================#
# Functions for running the Schisto sim!
#=================================================================================#

  # If larvae have been alive in human for > 35 days, turn into mating m/f adult worm
  function larvae_maturity!(human::Human,p::MyParameters)

      # if there are larvea that have been inside a host for >35 days,
      # mature into male/female adult worm
      if length(human.larvae) > 35

          # not used
          #bin = Binomial(human.larvae[1], 0.5 + 0.5*(human.vac_status>0)*0.975)

          bin = Binomial(human.larvae[1], 0.5)
          x = rand(bin)
          human.male_worms[1] = human.male_worms[1] + x
          human.female_worms[1] = human.female_worms[1] + human.larvae[1] - x

          # age the larvae population by a day
          shift!(human.larvae)
      end
  end

  # Adult worms can be in a number 'worm_stages' stages of life lasting on average 1/prob days
  # This is done to generate gamma distributed lifespans
  function worm_maturity!(human::Human,p::MyParameters)

      # hard coded in worm lifespan, to be changed later
      prob = p.worm_stages / (5.7*365)

      # if worms are in last stage of life, a certain fraction of them die
      bin_m = Binomial(human.male_worms[p.worm_stages], 1.0-prob)
      bin_f = Binomial(human.female_worms[p.worm_stages], 1.0-prob)
      human.male_worms[p.worm_stages] = rand(bin_m)
      human.female_worms[p.worm_stages] = rand(bin_f)

      # cycling back through worm stages, so that worms cant travel through the whole array in 1 day
      for i in (p.worm_stages-1):-1:1
          bin_m = Binomial(human.male_worms[i], prob)
          bin_f = Binomial(human.female_worms[i], prob)
          aging_males = rand(bin_m)
          aging_females = rand(bin_f)

          # test to make sure sampling correctly
          if typeof(aging_males) == Float64
              print(aging_males)
          end

          if typeof(aging_females) == Float64
              print(aging_females)
          end

          human.male_worms[i+1] = human.male_worms[i+1] + aging_males
          human.male_worms[i] = human.male_worms[i] - aging_males
          human.female_worms[i+1] = human.female_worms[i+1] + aging_females
          human.female_worms[i] = human.female_worms[i] - aging_females
      end


  end

  # total number of possible worm pairs in a host
  function worm_pairs(human::Human)
    return min(sum(human.male_worms), sum(human.female_worms))
  end

  # egg production
  function eggs_produced!(human::Human,p::MyParameters,v::Vaccination_Program)

      wp = worm_pairs(human)

      # alternative fit to cheever egg data
      '''
      mean = 100.996 * wp^0.458127 / (1.0 + 2.46988* exp(-0.142443*wp))
      eggs  = Poisson((1.0  - 1.0(human.vac_status>0))*mean/3.0)
      human.eggs = round(rand(eggs))
      '''
      human.eggs = round((1.0 - (human.vac_status>0)*v.fecundity_eff) * p.max_fecundity * wp * exp(-0.0007 * (sum(human.male_worms)+sum(human.female_worms))))

  end

  # loads death rates per age from file
  function load_age_death_rates()

    # use the default demographic data
    demo_file = open(readdlm,"data/Demographies2.txt")

    L = 0
    for i=1:length(demo_file[2,:])
      if (typeof(demo_file[2,i]) == Int64)
        L = L + 1
      end
    end
    age_upper_bounds = Array{Int64}(demo_file[2,2:L])
    deaths_per_year = Array{Float64}(demo_file[1,2:L])

    # transform deaths per year to deaths per day
    deaths_per_day = zeros(deaths_per_year)
    deaths_per_day = 1.0 - (1.0.-deaths_per_year).^(1.0/365.0)
    #deaths_per_day = 1.0 - exp.(-1.0*deaths_per_year/365.0)
    human_death_rates = hcat(age_upper_bounds,deaths_per_day)

    return human_death_rates
  end

  # load contact_rates from a given file
  function load_contact_data(file::String)
    return open(readdlm,file)
  end

  function death_rate(age::Int,p::MyParameters)

    for j in 1:length(p.death_human[:,1])
      if age < 365*p.death_human[j,1]  # convert from years to days
        return p.death_human[j,2]
      end
    end

    return p.death_human[length(p.death_human[:,1]),2]
  end


  # kill humans with given rates, and make sure children are born to make up the numbers
  function human_deaths!(e::Environment,p::MyParameters)

    # Systematic non compliance beta-params
    a = p.snc_coverage * (1.0 - p.snc_correlation) / p.snc_correlation
    b = (1.0 - p.snc_coverage)*(1.0-p.snc_correlation) / p.snc_correlation
    beta = Distributions.Beta(a,b)

    new_borns = Array{Human,1}(0)
    n_humans_pre = length(e.humans)

    # For every person, working backwards so we can iterate correctly
    for i in n_humans_pre:-1:1

      # if human older than biggest age in file, try to kill them with high prob
      if e.humans[i].age > 365*p.death_human[length(p.death_human[:,1]),1]
        if rand() < 50.0/(1*365.0)
          splice!(e.humans,i) # kill!

          new_human = Human() # human birth
          new_human.age = 1
          new_human.male_worms = zeros(p.worm_stages)         # change later to account for different worm_stage nos
          new_human.female_worms = zeros(p.worm_stages)         #
          new_human.predisposition = rand(p.gamma_r)
          new_human.snc_prob = rand(beta)
          push!(new_borns,new_human)
        end
        continue
      end

      # find age in death rates
      for j in 1:length(p.death_human[:,1])

        if e.humans[i].age < 365*p.death_human[j,1]  # 365 is to convert from years to days
          if rand() < p.death_human[j,2]
            splice!(e.humans,i) # kill!

            new_human = Human()
            new_human.age = 1
            new_human.male_worms = zeros(p.worm_stages)          # change later to account for different worm_stage nos
            new_human.female_worms = zeros(p.worm_stages)
            new_human.predisposition = rand(p.gamma_r)
            new_human.snc_prob = rand(beta)
            push!(new_borns,new_human)
          end
          break
        end
      end
    end

    # add newborns to survived population
    e.humans = vcat(e.humans,new_borns)

    # if something has wrong...
    n_humans_post = length(e.humans)
    difference = n_humans_pre - n_humans_post
    if length(e.humans)!=n_humans_pre
      println("shit")
    end

  end

  # turn eggs into larvae in environment
  function larvae_production!(e::Environment,p::MyParameters)
      new_larvae = 0
      for human in e.humans
          new_larvae += round(1.0*human.eggs*age_dependent_contact(div(human.age,365),p.age_contact))
      end
      new_larvae = round(new_larvae)
      if new_larvae < 0
        println("uh oh")
      end
      push!(e.larvae,new_larvae)
  end

  # find age dependent contact rate for a given age
  function age_dependent_contact(age::Int64,age_contacts::Array{Float64,2})
    L =  size(age_contacts)[1]
    for i in 1:L
      if age<age_contacts[i,1]
        return age_contacts[i,2]
      end
    end
    return age_contacts[L,2]
  end

  # human uptake from environment
  function larvae_uptake!(e::Environment,p::MyParameters,shuffle_pop::Int64,v::Vaccination_Program)

    # if larvae have lived for 40 days, uptake
    if length(e.larvae) > 40

      # should shuffle pop just to be safe (when transmission is low), but maybe this slows it down?
      if shuffle_pop == 1
        shuffle!(e.humans)
      end

      # humans get infected by larvae in their last two days of life
#      for human in e.humans

#        mean1 = human.predisposition*p.contact_rate*age_dependent_contact(div(human.age,365),p.age_contact)*e.larvae[1]
#        mean2 = human.predisposition*p.contact_rate*age_dependent_contact(div(human.age,365),p.age_contact)*e.larvae[2]
#        mean3 = human.predisposition*p.contact_rate*age_dependent_contact(div(human.age,365),p.age_contact)*e.larvae[3]
#        mean1 = mean1*(1 - (human.vac_status>0)*v.contact_eff)
#        mean2 = mean2*(1 - (human.vac_status>0)*v.contact_eff)
#        mean3 = mean3*(1 - (human.vac_status>0)*v.contact_eff)
#
#        uptake1 = rand(Poisson(mean1))
#        uptake2 = rand(Poisson(mean2))
#        uptake3 = rand(Poisson(mean3))
#        e.larvae[1] = e.larvae[1] - uptake1
#        e.larvae[2] = e.larvae[2] - uptake2
#        e.larvae[3] = e.larvae[3] - uptake3

        # for vaccine
        #push!(human.larvae,round((l1+l2)*(1.0-human.vacc_protect)))
        #push!(human.larvae,round((uptake1+uptake2)*(1.0-human.vacc_protect*0.5)))

#        push!(human.larvae,round((uptake1+uptake2+uptake3)))
#        end
      if e.larvae[1] < 0      # error check
        println(e.larvae[1])
      end
      e.infective_larvae = e.infective_larvae + shift!(e.larvae)
    end

    for human in e.humans

      mean1 = human.predisposition*p.contact_rate*age_dependent_contact(div(human.age,365),p.age_contact)*e.infective_larvae
      mean1 = mean1*(1.0 - (human.vac_status>0)*v.contact_eff) # uptake dependent on vac status

      if mean1 < 0 # error check
        println(human.predisposition)
        println(p.contact_rate)
        println(age_dependent_contact(div(human.age,365),p.age_contact))
        println(e.infective_larvae)
      end

      uptake1 = rand(Poisson(mean1))
      push!(human.larvae,round(uptake1))
      e.infective_larvae = max(0,e.infective_larvae - uptake1)
    end

    e.infective_larvae = Int(max(0,round(e.infective_larvae*(1.0 - 1.0/1.4))))  # infective worms dying
  end


  # vaccine runs out after x days
  function vac_decay!(human::Human)
    if human.vac_status > 0
      human.vac_status = human.vac_status - 1
    end
  end


  # do all the things needed to do in a normal day
  function update_env!(e::Environment,p::MyParameters,v::Vaccination_Program)
      e.time += 1
      for human in e.humans
          human.age = human.age + 1
          larvae_maturity!(human,p)
          eggs_produced!(human,p,v)
          worm_maturity!(human,p)
          vac_decay!(human)
      end
      larvae_production!(e,p)
      human_deaths!(e,p)
      larvae_uptake!(e,p,1,v)
  end


  # create population of hosts
  function create_population(p::MyParameters,N::Int64,larvae0::Int64,worms0::Int64)
      env = Environment()
      u = Geometric(1.0/(20*365))
      a = p.snc_coverage * (1.0 - p.snc_correlation) / p.snc_correlation
      b = (1.0 - p.snc_coverage)*(1.0-p.snc_correlation) / p.snc_correlation
      beta = Distributions.Beta(a,b)

      for i in 1:N
          new_human = Human()
          new_human.age = round(1+rand(u))
          new_human.male_worms = zeros(p.worm_stages)
          new_human.male_worms[1] = worms0
          new_human.female_worms = zeros(p.worm_stages)
          new_human.female_worms[1] = worms0
          new_human.predisposition = rand(p.gamma_r)
          new_human.snc_prob = rand(beta)
          push!(env.humans,new_human)
      end

      for i in 1:30
          push!(env.larvae,larvae0)
      end

      return env
  end

  # histograms collected from data, non-essential
  function get_distributions(env::Environment)

      eggs = zeros(length(env.humans))
      burden = zeros(length(env.humans))
      age = zeros(length(env.humans))
      prevalence = zeros(length(env.humans))

      for i in 1:length(env.humans)
          burden[i] = 1.0*sum(env.humans[i].male_worms) + sum(env.humans[i].female_worms)
          age[i] = 1.0*env.humans[i].age
          eggs[i] = 1.0*env.humans[i].eggs
      end

      df = DataFrame(Age = 1:length(env.humans),Burden = 1:length(env.humans), Eggs = 1:length(env.humans))
      bins = collect(0.0:3.0:2*99.0)
      ages = cut(round.(1.0*age/365.0),bins)
      df[:Age] = ages;
      df[:Burden] = burden;
      df[:Eggs] = eggs;

      burden_histogram = by(df, [:Age], df -> mean(df[:Burden]));
      egg_histogram = by(df, [:Age], df -> mean(df[:Eggs]));
      age_histogram = by(df, [:Age], df -> length(df[:Eggs]));
      prevalence_histogram = by(df, [:Age], df -> sum(df[:Eggs].>0)/length(df[:Burden]))

      dist = Epi_distributions(age_histogram,burden_histogram,egg_histogram,prevalence_histogram,1)
      return dist
  end

  function get_distributions_byyear(env::Environment)

      eggs = zeros(length(env.humans))
      burden = zeros(length(env.humans))
      age = zeros(length(env.humans))
      prevalence = zeros(length(env.humans))

      for i in 1:length(env.humans)
          burden[i] = 1.0*sum(env.humans[i].male_worms) + sum(env.humans[i].female_worms)
          age[i] = 1.0*env.humans[i].age
          eggs[i] = 1.0*env.humans[i].eggs
      end

      df = DataFrame(Age = 1:length(env.humans),Burden = 1:length(env.humans), Eggs = 1:length(env.humans))
      bins = collect(0.0:1.0:2*140.0)
      ages = cut(round.(1.0*age/365.0),bins)
      df[:Age] = ages;
      df[:Burden] = burden;
      df[:Eggs] = eggs;

      burden_histogram = by(df, [:Age], df -> mean(df[:Burden]));
      egg_histogram = by(df, [:Age], df -> mean(df[:Eggs]));
      age_histogram = by(df, [:Age], df -> length(df[:Eggs]));
      prevalence_histogram = by(df, [:Age], df -> sum(df[:Eggs].>0)/length(df[:Burden]))

      dist = Epi_distributions(age_histogram,burden_histogram,egg_histogram,prevalence_histogram,1)
      return dist
  end

  # update histograms
  function update_distributions!(dists::Epi_distributions, env::Environment)
      new_dists = get_distributions(env)

      age_cat = vcat(new_dists.age_histogram,dists.age_histogram)
      burden_cat = vcat(new_dists.burden_histogram,dists.burden_histogram)
      egg_cat = vcat(new_dists.egg_histogram,dists.egg_histogram)
      prev_cat = vcat(new_dists.prevalence_histogram,dists.prevalence_histogram)

      dists.age_histogram = by(age_cat,[:Age], df->sum(df[:x1]) )
      dists.burden_histogram = by(burden_cat,[:Age], df->sum(df[:x1]) )
      dists.egg_histogram = by(egg_cat,[:Age], df->sum(df[:x1]) )
      dists.prevalence_histogram = by(prev_cat,[:Age], df->sum(df[:x1]) )

      dists.n_samples = dists.n_samples + 1
  end

  function update_distributions_byyear!(dists::Epi_distributions, env::Environment)
      new_dists = get_distributions_byyear(env)

      age_cat = vcat(new_dists.age_histogram,dists.age_histogram)
      burden_cat = vcat(new_dists.burden_histogram,dists.burden_histogram)
      egg_cat = vcat(new_dists.egg_histogram,dists.egg_histogram)
      prev_cat = vcat(new_dists.prevalence_histogram,dists.prevalence_histogram)

      dists.age_histogram = by(age_cat,[:Age], df->sum(df[:x1]) )
      dists.burden_histogram = by(burden_cat,[:Age], df->sum(df[:x1]) )
      dists.egg_histogram = by(egg_cat,[:Age], df->sum(df[:x1]) )
      dists.prevalence_histogram = by(prev_cat,[:Age], df->sum(df[:x1]) )

      dists.n_samples = dists.n_samples + 1
  end

  # rough calculation of r0
  function calculateR0(p::MyParameters,env::Environment,dists::Epi_distributions)

    n = length(env.humans)
    #age = dists.age_histogram[:x1][1:23]/dists.n_samples
    age = dists.age_histogram[:x1]/dists.n_samples
    #println(dists.age_histogram)

    bins = collect(3.0:3.0:(length(age))*3.0)
    con = [age_dependent_contact(Int(i),p.age_contact) for i in bins]
    g = Distributions.Gamma(p.worm_stages,(5.7/p.worm_stages)*365.0)


    d_r = [1.0-death_rate(Int(365*bins[i]),p) for i in 1:length(bins)]
    d_r = -365.0*log.(d_r)
    d_r = [death_rate(Int(365*bins[i]),p) for i in 1:length(bins)]

    integrand = [ 3.0 * d_r[i] * age[i] * sum([3.0 * con[j]*(1.0 - cdf(g,365*(bins[i] - bins[j]))) for j in 1:i]) for i in 1:length(bins)]
    d_h = 1.0/(365.0*20.0)
    d_h = sum([d_r[i]*age[i] for i in 1:length(bins)])/n
    #d_h = mean(d_r)
    #println(1.0/(d_h*365.0))
    d_w = 1.0/(365.0*5.7)
    d = 1.0/1.40
    lamb = 1.0*p.max_fecundity*sum(con.*age)/n
    b = p.contact_rate*sum(con.*age)/n
    b2 = p.contact_rate * sum(integrand)
    #println("psi is: ",   (0.5*p.max_fecundity*exp(-0.0007))/((d_h+d_w)*d))

    #return 0.5*b*n*p.max_fecundity/((d_h+d_w)*(d+n*p.contact_rate))
    #return 0.5*b*n*p.max_fecundity*exp(-0.0007)/((d_h+d_w)*(d))
    #return 1.0 *(2.0/3.0)* b2 * n * p.max_fecundity*exp(-0.0007) / d
    0.5*b*n*lamb*exp(-0.0007)/((d_h+d_w)*(d))
  end

  function set_beta(R0::Float64,p::MyParameters,env::Environment,dists::Epi_distributions)

    n = length(env.humans)
    l = length(dists.age_histogram[:x1])
    age = dists.age_histogram[:x1][1:l]/dists.n_samples
    bins = collect(3.0:3.0:(length(age))*3.0)
    con = [age_dependent_contact(Int(i),p.age_contact) for i in bins]
    g = Distributions.Gamma(p.worm_stages,(5.7/p.worm_stages)*365.0)
    d_r = [1.0-death_rate(Int(365*bins[i]),p) for i in 1:length(bins)]
    d_r = -365.0*log.(d_r)
    d_r = [death_rate(Int(365*bins[i]),p) for i in 1:length(bins)]

    integrand = [ 3.0 * d_r[i] * age[i] * sum([3.0 * con[j]*(1.0 - cdf(g,365*(bins[i] - bins[j]))) for j in 1:i]) for i in 1:length(bins)]

    d_h = 1.0/(365.0*20.0)
    d_h = sum([d_r[i]*age[i] for i in 1:length(bins)])/n
    d_w = 1.0/(5.7*365.0)
    d = 1.0/1.4
    lamb = 1.0*p.max_fecundity*sum(con.*age)/n

    #b = R0*(d_h+d_w)*d / (0.5*n*p.max_fecundity*exp(-0.0007))
    b = R0*(d_h+d_w)*d / (0.5*n*lamb*exp(-0.0007))

    b2 = R0*d / (1.0n*p.max_fecundity*exp(-0.0007))
    return b*n/sum(con.*age)
    #return b2 / sum(integrand)
  end


  # calculate aggregation statistic from population
  function aggregation(e::Environment)
      burden = zeros(length(e.humans))
      for i in 1:length(burden)
          burden[i] = sum(e.humans[i].female_worms) + sum(e.humans[i].male_worms)
      end

      if sum(burden) > 1
          return mean(burden)^2 / (var(burden) - mean(burden))
      else
          return 0.0
      end
  end

  # calculate mean_age statistic from pop
  function mean_age(e::Environment)
      a=0.0
      for human in e.humans
          a = a + human.age
      end
      return a/(length(e.humans)*365.0)
  end

  # simulate the kato katz test, average from two "slides"
  function kato_katz(p::MyParameters,h::Human)
      gamma1 = rand(p.gamma_k)
      gamma2 = rand(p.gamma_k)
      gamma3 = rand(p.gamma_k)
      pois1 = Poisson(gamma1*h.eggs)
      pois2 = Poisson(gamma2*h.eggs)
      pois3 = Poisson(gamma3*h.eggs)
      #return Int(floor(1.0(rand(pois1)))) #+rand(pois2))))   # Can use an average of two or three
      return (rand(pois1)+rand(pois2))/2.0
  end

  # calculate stratified prevalences, within age range
  function stratified_prevalence(age_low::Int64,age_high::Int64,e::Environment,p::MyParameters)

    high = 0
    medium = 0
    low = 0
    total = 0

    for human in e.humans
      if rand() < 1.5 # this is trivially always true, not sure why this is here
        if (human.age > age_low*365) & (human.age < age_high*365)

          total = total + 1
          kk_test = kato_katz(p,human)

          if kk_test >= 17
              high = high + 1
          elseif kk_test > 4
              medium = medium + 1
          elseif kk_test > 0
              low = low + 1
          end

        end
      end
    end

    if total > 0
      return 1.0*high/total,1.0*medium/total,1.0*low/total
    else
      return 0.0,0.0,0.0
    end
  end

  # create arrays for stats
  function initialise_statistics(n_samples::Int64,ensembles::Int64)

    mean_burden = zeros(n_samples,ensembles)
    prevalence = zeros(n_samples,ensembles)
    reservoir = zeros(n_samples,ensembles)
    aggregation = zeros(n_samples,ensembles)
    stratified_prevalences = zeros(n_samples,3,ensembles)
    adult_prevalence = zeros(n_samples,ensembles)
    adult_stratified_prevalences = zeros(n_samples,3,ensembles)
    return Statistics(mean_burden,prevalence,aggregation,reservoir,stratified_prevalences,adult_prevalence,adult_stratified_prevalences)

  end


  function sample_statistics!(stats::Statistics,env::Environment,t::Int64,ensemble::Int64,p::MyParameters)

    N = length(env.humans)

    # calculate mean burden
    b = 0
    for j in 1:N
        b += sum(env.humans[j].male_worms) + sum(env.humans[j].female_worms)
    end
    stats.mean_burden[t,ensemble] = b/N

    # calculate stratefied prevalences
    high, med, low = stratified_prevalence(5,15,env,p)
    stats.stratified_prevalence[t,1,ensemble] = low
    stats.stratified_prevalence[t,2,ensemble] = med
    stats.stratified_prevalence[t,3,ensemble] = high
    stats.prevalence[t,ensemble] = low + med + high

    ahigh, amed, alow = stratified_prevalence(15,140,env,p)
    stats.adult_stratified_prevalence[t,1,ensemble] = alow
    stats.adult_stratified_prevalence[t,2,ensemble] = amed
    stats.adult_stratified_prevalence[t,3,ensemble] = ahigh
    stats.adult_prevalence[t,ensemble] = alow + amed + ahigh

    stats.worms_in_pool[t,ensemble] = sum(env.larvae)
    stats.aggregation[t,ensemble] = aggregation(env)

  end

  function save_environment(name::String,e::Environment)
      l = length(e.humans)
      filename = name * string(l) * ".jld"
      save(filename, "environment", e)
  end

  function load_environment(name::String)
      d = load(name)
      return d["environment"]
  end

  function reset_snc!(env::Environment, p::MyParameters, cor::Float64, cov::Float64)

    if cov == 0.0
      for human in env.humans
        human.snc_prob = 0
      end
    else

      a = cov * (1.0 - cor) / cor
      b = (1.0 - cov)*(1.0-cor) / cor
      beta = Distributions.Beta(a,b)
      for human in env.humans
        human.snc_prob = rand(beta)
      end
    end
  end

  # apply MDA to a population
  function MDA!(env::Environment,mda_pro::MDA_Program,p::MyParameters)

      for human in env.humans
          ages = collect(0:1:200)
          if div(human.age,365) in [5,6,7,8,9,10,11,12,13,14]    # MDA aimed at SACs
          #if div(human.age,365) in ages
            if rand() < human.snc_prob
              for j in 1:p.worm_stages
                human.female_worms[j] = round((1.0 - mda_pro.efficacy)*human.female_worms[j])
                human.male_worms[j] = round((1.0 - mda_pro.efficacy)*human.male_worms[j])
              end
              human.eggs = 0      # no eggs shed today
              human.larvae = []   # probably not neccessary
            end
          end
      end
  end

# apply vaccination to population
  function vaccinate!(env::Environment,v_pro::Vaccination_Program,p::MyParameters)

    for human in env.humans
      if length(find(v_pro.target_ages .== div(human.age,365))) > 0    # is this person in list od ages to be vaccinated
        if rand() < v_pro.coverages[find(v_pro.target_ages .== div(human.age,365))][1]
          human.vac_status = Int(round(365*(v_pro.duration + v_pro.duration*randn()/10.0))) # x years of vaccination protection
          for j in 1:p.worm_stages
            human.female_worms[j] = round((1.0 - v_pro.theraputic_eff)*human.female_worms[j])
            human.male_worms[j] = round((1.0 - v_pro.theraputic_eff)*human.male_worms[j])
          end
          #human.larvae = []
          #human.eggs= 0
        end
      end
    end
  end

  # make some plots of the stats through time
  function make_plots(stats::Statistics,T::Int64,sample_rate::Int64)

    PyPlot.plt[:style][:use]("ggplot")          # chosen style for plotting
    x = (1.0/365.0)*linspace(1,T,T/sample_rate)
    fig = plt[:figure]()

    subplot(2,2,1)
    plot(x,stats.mean_burden,"grey",linewidth=0.2,alpha=0.4)
    plot(x,mean(stats.mean_burden,2),linewidth=0.7)
    ax = gca() # get current axes
    ax[:grid](linestyle="--", linewidth=0.5)
    title("Mean burden")

    #=
    subplot(2,2,4)
    plot(x,stats.aggregation,"grey",linewidth=0.2,alpha=0.4)
    plot(x,mean(stats.aggregation,2),linewidth=0.7)
    ax = gca() # get current axes
    ax[:grid](linestyle="--", linewidth=0.5)
    title("Aggregation")
    ax[:set_ylim]((0,0.6));
    xlabel("Years")
    =#

    subplot(2,2,4)
    plot(x,reshape(mean(stats.stratified_prevalence[:,3,:],2),(length(x))))
    ax = gca() # get current axes
    ax[:grid](linestyle="--", linewidth=0.5)
    title("Low/med/high burden Prevalence")
    xlabel("Years")
    ax[:set_ylim]((0.0,1.0));

    subplot(2,2,3)
    stackplot(x,reshape(mean(stats.stratified_prevalence[:,1,:],2),(length(x))),
        reshape(mean(stats.stratified_prevalence[:,2,:],2),(length(x))),
        reshape(mean(stats.stratified_prevalence[:,3,:],2),(length(x))),
        colors=["gold","darkorange","tomato"])
    ax = gca() # get current axes
    ax[:grid](linestyle="--", linewidth=0.5)
    title("Low/med/high burden Prevalence")
    xlabel("Years")
    ax[:set_ylim]((0.0,1.0));

    subplot(2,2,2)
    #plot(x,stats.worms_in_pool,"grey",linewidth=0.2,alpha=0.4)
    stackplot(x,reshape(mean(stats.adult_stratified_prevalence[:,1,:],2),(length(x))),
    reshape(mean(stats.adult_stratified_prevalence[:,2,:],2),(length(x))),
    reshape(mean(stats.adult_stratified_prevalence[:,3,:],2),(length(x))),
    colors=["gold","darkorange","tomato"])
    ax = gca() # get current axes
    ax[:grid](linestyle="--", linewidth=0.5)
    title("adult")
    xlabel("Years")
    ax[:set_ylim]((0.0,1.0));

    fig[:suptitle]("Schisto time plots")
    fig[:tight_layout]()
    fig[:subplots_adjust](top=0.88)
    savefig("equilibrium.png", bbox_inches="tight") # pick a better file name to save to
  end

  function make_prevalence_plot(stats::Statistics,T::Int64,sample_rate::Int64)

    PyPlot.plt[:style][:use]("ggplot")          # chosen style for plotting
    x = (1.0/365.0)*linspace(1,T,T/sample_rate) - 1.0
    fig = plt[:figure]()

    plot(x,stats.prevalence,"grey",linewidth=0.2,alpha=0.25)
    plot(x,mean(stats.prevalence,2),linewidth=0.7)
    ax = gca() # get current axes
    ax[:grid](linestyle="--", linewidth=0.5)
    xlabel("Years")
    ylabel("Prevalence")
    title("Prevalence")
    savefig("prevalence_plot.png", bbox_inches="tight") # pick a better file name to save to
  end

  function make_high_prevalence_plot(stats::Statistics,T::Int64,sample_rate::Int64)

    PyPlot.plt[:style][:use]("seaborn-darkgrid")          # chosen style for plotting
    x = (1.0/365.0)*linspace(1,T,T/sample_rate) - 1.0
    fig = plt[:figure]()
    y1 = 5*ones(x)
    y2 = 1*ones(x)

      #plot(x,stats.stratified_prevalence[:,3,:],"grey",linewidth=0.1,alpha=0.2)


      #plot((10,10),(0,0.2),linestyle="--",linewidth=0.85,color="steelblue")
      #plot((5,5),(0,0.2),linestyle="--",linewidth=0.85,color="steelblue")

      ax = gca() # get current axes

      # high burden SAC
      plot(x,100.0*mean(stats.stratified_prevalence[:,3,:],2),"orangered",linewidth=1.0)
      lb = 100.0*max.((mean(stats.stratified_prevalence[:,3,:],2) - 2.0*std(stats.stratified_prevalence[:,3,:],2))[:],0.0)
      ub = 100.0(mean(stats.stratified_prevalence[:,3,:],2) + 2.0*std(stats.stratified_prevalence[:,3,:],2))[:]
      ax[:fill_between](x,lb,ub,alpha=0.2,color="orangered")

      # high burden Adults
      plot(x,100.0*mean(stats.adult_stratified_prevalence[:,3,:],2),"royalblue",linewidth=1.0)
      lb = 100.0*max.((mean(stats.adult_stratified_prevalence[:,3,:],2) - 2.0*std(stats.adult_stratified_prevalence[:,3,:],2))[:],0.0)
      ub = 100.0*(mean(stats.adult_stratified_prevalence[:,3,:],2) + 2.0*std(stats.adult_stratified_prevalence[:,3,:],2))[:]
      ax[:fill_between](x,lb,ub,alpha=0.2,color="royalblue")

      # all SACS
      plot(x,100.0*mean(stats.prevalence,2),linewidth=1.0,linestyle="--",color="orangered")
      lb = 100.0*max.((mean(stats.prevalence,2) - 2.0*std(stats.prevalence,2))[:],0.0)
      ub = 100.0*(mean(stats.prevalence,2) + 2.0*std(stats.prevalence,2))[:]
      ax[:fill_between](x,lb,ub,alpha=0.2,color="orangered")

      # all adults
      plot(x,100.0*mean(stats.adult_prevalence,2),linewidth=1.0,linestyle="--",color="royalblue")
      lb = 100.0*max.((mean(stats.adult_prevalence,2) - 2.0*std(stats.adult_prevalence,2))[:],0.0)
      ub = 100.0*(mean(stats.adult_prevalence,2) + 2.0*std(stats.adult_prevalence,2))[:]
      ax[:fill_between](x,lb,ub,alpha=0.2,color="royalblue")

      plot(x,y1,linewidth=0.85,linestyle="--",color="steelblue")
      plot(x,y2,linewidth=0.85,linestyle="--",color="steelblue")


      ax[:set_yticks](linspace(0,100,11), minor=true)
      ax[:set_xticks](linspace(0,15,16), minor=true)
      ax[:grid](b=true, which="minor",color="w", linewidth=0.5,linestyle="--")
      upper_lim = 100.0 * min(1.04,1.33 * mean(stats.prevalence,2)[1])
      println(upper_lim)
      ax[:set_ylim]((-0.4,upper_lim));
      ax[:set_xlim]((-0.5,div(T,365) - 1.0 + 0.2));
      ax[:grid](linestyle="--", linewidth=0.5)
      xlabel("Time (Years)")
      ylabel("Prevalence (%)")

  #ax[:get_xaxis][:set_minor_locator](mpl.ticker.AutoMinorLocator())
  #ax[:get_yaxis][:set_minor_locator](mpl.ticker.AutoMinorLocator())
  #ax[:grid](b=true, which="major", color="w", linewidth=1.0)
  #ax[:grid](b=true, which="minor")


  #ax[:set_yticks]([-1.25, -0.75, -0.25,0.24,0.75,1.25], minor=true)

  #plt[:figlegend](["a","a","a","a"],loc="upper right")

  box = ax[:get_position]()
  #ax[:set_position]([box[:x0], box[:y0] + box[:height] * 0.1,
  #                box[:width], box[:height] * 0.9])
  ax[:legend](["SAC high-burden prevalence","Adult high-burden prevalence","SAC prevalence","Adult prevalence"],loc="upper right",# bbox_to_anchor=(0.5, -0.1),
            fancybox=true, shadow=true, ncol=1)
  title("Ages: 5, 10 years; Coverage 75%, 60%; Duration: 5 years")
    savefig("high_prevalence_plot.png", bbox_inches="tight") # pick a better file name to save to
  end



  function make_prevalence_histogram(dists1::Epi_distributions,dists2::Epi_distributions)
    x = 3:3:25*3
    PyPlot.plt[:style][:use]("seaborn-darkgrid")
    fig1 = plt[:figure]()
    PyPlot.bar(dists1.egg_histogram[:Age][1:25],24.0*dists1.egg_histogram[:x1][1:25]/dists1.n_samples,alpha=0.5,label="Pre")
    PyPlot.bar(dists2.egg_histogram[:Age][1:25],24.0*dists2.egg_histogram[:x1][1:25]/dists2.n_samples,alpha=0.5,label="Post",color="green")
    PyPlot.legend(loc="upper right")
    ax = gca()
    ax[:set_xticklabels](x)
    title("Egg Output")
    xlabel("Age")
    ylabel("Eggs per gram")
  end

  # run a certain number of time
  function run_timesteps!(ensembles::Int64, T::Int64, sample_rate::Int64,env0::Environment,stats::Statistics,p::MyParameters,v::Vaccination_Program)

    # create environment
    env = Environment()
    pre_dists = Array{Epi_distributions,1}(0)
    post_dists = Array{Epi_distributions,1}(0)

    # bring up progress bar, not really neccessary (@ symbol means macro)
    @showprogress 1 "Progress..." for k in 1:ensembles
        #srand(k*100)
        srand()

        # Make a copy of initial condition, so that we can take multiple run from same initial condition
        env = deepcopy(env0)
        env.time = 1
        distributions = get_distributions(env)
        initial_distributions = get_distributions(env)
        final_distributions = get_distributions(env)

        # run through days
        for i in 1:T

            update_env!(env,p,v)

            # sample statistics with given rate
            if i%sample_rate == 0
                sample_statistics!(stats,env,div(i,sample_rate),k,p)
                if stats.prevalence[div(i,sample_rate),k] < 0.01
                  break
                end
            end

            if (i%20 == 0) & (i< T/20)
              update_distributions!(initial_distributions,env)
            end

            if (i%20 == 0) & (i> 19*T/20)
              update_distributions!(final_distributions,env)
            end

            # update distributions every year
            if (i%365 == 0) & (i> T/2)
              update_distributions!(distributions,env)
            end

        end

        r0 = calculateR0(p,env,distributions)
        push!(pre_dists,initial_distributions)
        push!(post_dists,final_distributions)
        #make_prevalence_histogram(distributions)
        #println("final R0 is ", r0)
    end

    return env, pre_dists,post_dists
  end

  # vaccination routine
  function run_timesteps_vaccine!(ensembles::Int64, T::Int64, sample_rate::Int64,env0::Environment,stats::Statistics,p::MyParameters,v_pro::Vaccination_Program)

    pre_dists = Array{Epi_distributions,1}(0)
    post_dists = Array{Epi_distributions,1}(0)
    frequency = 1  # how many years between each vac, shouldnt really be hardcoded

    @showprogress 1 "Progress..." for k in 1:ensembles
        srand()

        # Comment out as desired
        env = deepcopy(env0)
        L = rand(1:365*16)
        v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")

        for i in 1:L
            update_env!(env,p,v)
        end

        env.time = 1
        distributions = get_distributions(env)
        initial_distributions = get_distributions(env)
        final_distributions = get_distributions(env)

        for i in 1:T
            update_env!(env,p,v_pro)
            if i%sample_rate == 0
                sample_statistics!(stats,env,div(i,sample_rate),k,p)
            end

            #if (i%365 == 0) #& (i < 365*10)
            if (( ((i-365)/365)%frequency) == 0)
              vaccinate!(env,v_pro,p)
            end

            if (i%5 == 0) & (i< 1*T/16)
              update_distributions!(initial_distributions,env)
            end

            if (i%5 == 0) & (i> 15*T/16)
              update_distributions!(final_distributions,env)
            end
        end
        push!(pre_dists,initial_distributions)
        push!(post_dists,final_distributions)
    end

    return stats, pre_dists, post_dists
  end

  # vaccination with a catchup campaign at year 1
  function run_timesteps_vcatchup!(ensembles::Int64, T::Int64, sample_rate::Int64,env0::Environment,stats::Statistics,p::MyParameters,v_pro::Vaccination_Program)

    pre_dists = Array{Epi_distributions,1}(0)
    post_dists = Array{Epi_distributions,1}(0)

    a = collect(1.0:1.0:15.0)
    cov = coverages(a)
    catchup = Vaccination_Program(v_pro.duration,a,cov,16,1.0,1.0,0.0,0.0,v_pro.intensity)

    @showprogress 1 "Progress..." for k in 1:ensembles
        srand()

        # Comment out as desired
        env = deepcopy(env0)

        L = rand(1:365*16) # gives random burn in time to decorrelate ICs a bit
        v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")
        for i in 1:L
            update_env!(env,p,v)
        end
        env.time = 1
        distributions = get_distributions(env)
        initial_distributions = get_distributions(env)
        final_distributions = get_distributions(env)

        for i in 1:T
            update_env!(env,p,v_pro)
            if i%sample_rate == 0
                sample_statistics!(stats,env,div(i,sample_rate),k,p)
            end

            if (i%365 == 0) & (i > 365)
              vaccinate!(env,v_pro,p)
            end

            if (i == 365)
              vaccinate!(env,catchup,p)
            end

            if (i%5 == 0) & (i< 1*T/16)
              update_distributions!(initial_distributions,env)
            end

            if (i%5 == 0) & (i> 15*T/16)
              update_distributions!(final_distributions,env)
            end
        end
        push!(pre_dists,initial_distributions)
        push!(post_dists,final_distributions)
    end

    return stats, pre_dists, post_dists
  end

  # vaccination and MDA combined
  function run_timesteps_vacmda!(ensembles::Int64, T::Int64, sample_rate::Int64,env0::Environment,stats::Statistics,p::MyParameters,vacmda_pro::VACMDA_Program)

    pre_dists = Array{Epi_distributions,1}(0)
    post_dists = Array{Epi_distributions,1}(0)
    @showprogress 1 "Progress..." for k in 1:ensembles
        srand()

        # Comment out as desired
        env = deepcopy(env0)
        L = rand(1:365*16)
        v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")
        for i in 1:L
            update_env!(env,p,v)
        end
        env.time = 1
        distributions = get_distributions(env)
        initial_distributions = get_distributions(env)
        final_distributions = get_distributions(env)

        v_pro = Vaccination_Program(
          vacmda_pro.duration,
          vacmda_pro.target_ages,
          vacmda_pro.coverages,
          vacmda_pro.years,
          vacmda_pro.contact_eff,
          vacmda_pro.fecundity_eff,
          vacmda_pro.death_eff,
          vacmda_pro.theraputic_eff,
          vacmda_pro.intensity)

          ages = collect(1:1:150)
          cov = 1.0*zeros(ages)
          for i in 1:length(ages)
            if ages[i] < 5
              cov[i] = 0.0
            elseif ages[i] < 15
              cov[i] = 0.6
            else
              cov[i] = 0
            end
          end

        if vacmda_pro.intensity == "high"
          mda_pro = MDA_Program(cov,16,0.86,0.0,1)
        end

        if vacmda_pro.intensity == "medium"
          mda_pro = MDA_Program(cov,16,0.86,0.0,2)
        end

        if vacmda_pro.intensity == "low"
          mda_pro = MDA_Program(cov,16,0.86,0.0,3)
        end

        a = collect(1.0:1.0:15.0)
        cov = coverages(a)
        catchup = Vaccination_Program(v_pro.duration,a,cov,16,1.0,1.0,0.0,0.0,v_pro.intensity)


        for i in 1:T
            update_env!(env,p,v_pro)
            if i%sample_rate == 0
                sample_statistics!(stats,env,div(i,sample_rate),k,p)
            end

            if (( ((i-365)/365)%mda_pro.frequency) == 0) #& (i<365*11)
              reset_snc!(env, p, 0.000001, 0.6)
              MDA!(env,mda_pro,p)
            end

            #if (i%365 == 0)
            #  vaccinate!(env,v_pro,p)
            #end

            if (i%365 == 0) & (i > 365)
              vaccinate!(env,v_pro,p)
            end

            if (i == 365)
              vaccinate!(env,catchup,p)
            end

            if (i%5 == 0) & (i< 1*T/16)
              update_distributions!(initial_distributions,env)
            end

            if (i%5 == 0) & (i> 15*T/16)
              update_distributions!(final_distributions,env)
            end
        end
        push!(pre_dists,initial_distributions)
        push!(post_dists,final_distributions)
    end

    return stats, pre_dists, post_dists
  end



  # mda routine
  function run_timesteps_MDA!(ensembles::Int64, T::Int64, sample_rate::Int64,env0::Environment,stats::Statistics,p::MyParameters,mda_pro::MDA_Program,v::Vaccination_Program)

    pre_dists = Array{Epi_distributions,1}(0)
    post_dists = Array{Epi_distributions,1}(0)
    @showprogress 1 "Progress..." for k in 1:ensembles
        srand()

        # Comment out as desired
        env = deepcopy(env0)

        L = rand(1:365*16)
        vpro = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")
        for i in 1:L
            update_env!(env,p,vpro)
        end

        env.time = 1
        distributions = get_distributions(env)
        initial_distributions = get_distributions(env)
        final_distributions = get_distributions(env)
        #reset_snc!(env,p,0.25,0.6


        for i in 1:T
            update_env!(env,p,v)
            if i%sample_rate == 0
                sample_statistics!(stats,env,div(i,sample_rate),k,p)
            end
            if (( ((i-365)/365)%mda_pro.frequency) == 0) #& (i<365*11)
              reset_snc!(env, p, 0.4, mda_pro.coverage[10])
              MDA!(env,mda_pro,p)
            end

            if (i%10 == 0) & (i< 1*T/16)
              update_distributions!(initial_distributions,env)
            end

            if (i%10 == 0) & (i> 15*T/20)
              update_distributions!(final_distributions,env)
            end
        end

        push!(pre_dists,initial_distributions)
        push!(post_dists,final_distributions)
    end

    return stats, pre_dists, post_dists
  end

  # run to equilibrium
  function run_to_equilibrium!(env::Environment,n_hosts::Int64,ensembles::Int64,save::Bool,p::MyParameters,years::Int64)

    sample_rate = 25
    T = 365*years
    n_samples = div(T,sample_rate) # this is just integer division
    v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")

    p.death_human = load_age_death_rates()
    #p.age_contact = load_contact_data("data/IietuneContact.txt")
    #p.age_contact = load_contact_data("data/Iietune2.txt")
    p.age_contact = load_contact_data("data/ICL.txt")
    distributions = get_distributions(env)
    statistics = initialise_statistics(n_samples,ensembles)

    # run timesteps
    env, d1, d2 = run_timesteps!(ensembles,T,sample_rate,env,statistics,p,v)

    if save == true
      save_environment("environment", env)
    end

    r0 = calculateR0(p,env,d2[1])
    println("new R0 is ", r0)
    make_prevalence_histogram(d1[1],d2[1])

    return statistics, env, d2
  end

  function run_vaccination(env::Environment,v_pro::Vaccination_Program,p::MyParameters,dists::Epi_distributions)
    println("im working!")
    v_pro.years = 16
    years = v_pro.years
    sample_rate = years
    T = years*365
    ensembles = 250
    n_samples = div(T,sample_rate)

    statistics = initialise_statistics(n_samples,ensembles)

    if v_pro.intensity == "high"

      env = load_environment("ENV1000HIGH5.200e-4.jld")
      p.contact_rate = 5.200e-4
      p.gamma_r = Distributions.Gamma(0.24,1.0/0.24)
      p.max_fecundity = 0.14

    end

    if v_pro.intensity == "medium"
      env = load_environment("ENV1000MED1.543e-4.jld")
      p.contact_rate = 1.543e-4
      p.gamma_r = Distributions.Gamma(0.24,1.0/0.24)
      p.max_fecundity = 0.14

    end

    if v_pro.intensity == "low"

      env = load_environment("ENV1000LOW1.390e-4.jld")
      p.contact_rate = 1.390e-4
      p.gamma_r = Distributions.Gamma(0.04,1.0/0.04)
      p.max_fecundity = 0.14

    end

    v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"pro")
    for i in 1:365*1
        update_env!(env,p,v)
    end

    #env = load_environment("environment_HIGHPREV1000.jld")
    #p.contact_rate = set_beta(v_pro.R0,p,env,dists)

    catch_up = false

    if catch_up == true
      statistics, pre_dists,post_dists  = run_timesteps_vcatchup!(ensembles,T,sample_rate,env,statistics,p,v_pro)
    else
      statistics, pre_dists,post_dists  = run_timesteps_vaccine!(ensembles,T,sample_rate,env,statistics,p,v_pro)
    end

    pre = zeros(25)
    post = zeros(25)

    for i in 1:ensembles
      sorted_pre = sort(pre_dists[i].egg_histogram,cols = :Age)
      sorted_post = sort(post_dists[i].egg_histogram,cols = :Age)
      pre = pre + (1.0*convert(Array,sorted_pre[:x1][1:25]) / pre_dists[i].n_samples)/ensembles
      post = post + (1.0*convert(Array,sorted_post[:x1][1:25]) / post_dists[i].n_samples)/ensembles
    end

#    pre =  convert(Array,pre_dists[1].egg_histogram[:x1][1:25]) / pre_dists[1].n_samples
#    post = convert(Array,post_dists[1].egg_histogram[:x1][1:25]) / post_dists[1].n_samples#
#    println()
#    println(sort(pre_dists[1].egg_histogram,cols = :Age))
#    println(pre)
#    println()

    out = Output(statistics,pre,post)

    #=
    s = string(v_pro.target_ages)
    s = replace(s,", ",".")
    s = replace(s,"[","")
    s = replace(s,"]","")
    s = s*".duration."*string(Int(v_pro.duration))*"."*"ef"*string((v_pro.contact_eff^2))*".csv"

    ds = generate_age_dist(false);
    ages = convert(Array,ds.age_histogram[:x1]/ds.n_samples)
    get_economic_figs(statistics,v_pro.target_ages,v_pro.coverages,s,ages)

    #return convert(Array,post_dists[1].age_histogram[:,2])
    return 1.0
=#

    return out
  end

  function run_vacmda(env::Environment,vacmda_pro::VACMDA_Program,p::MyParameters,dists::Epi_distributions)
    println("im working!")
    years = vacmda_pro.years
    sample_rate = years
    T = years*365
    ensembles = 250
    n_samples = div(T,sample_rate)

    statistics = initialise_statistics(n_samples,ensembles)

    if vacmda_pro.intensity == "high"

      env = load_environment("ENV1000HIGH5.200e-4.jld")
      p.contact_rate = 5.00e-4
      p.gamma_r = Distributions.Gamma(0.24,1.0/0.24)
      p.max_fecundity = 0.14


    end

    if vacmda_pro.intensity == "medium"
      env = load_environment("ENV1000MED1.543e-4.jld")
      p.contact_rate = 1.543e-4
      p.gamma_r = Distributions.Gamma(0.24,1.0/0.24)
      p.max_fecundity = 0.14

    end

    if vacmda_pro.intensity == "low"
      env = load_environment("ENV1000LOW1.390e-4.jld")
      p.contact_rate = 1.390e-4
      p.gamma_r = Distributions.Gamma(0.04,1.0/0.04)
      p.max_fecundity = 0.14

    end


    v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")
    for i in 1:365*1
        update_env!(env,p,v)
    end

    #env = load_environment("environment_HIGHPREV1000.jld")
    #p.contact_rate = set_beta(v_pro.R0,p,env,dists)
    statistics, pre_dists,post_dists  = run_timesteps_vacmda!(ensembles,T,sample_rate,env,statistics,p,vacmda_pro)

    pre = zeros(25)
    post = zeros(25)

    for i in 1:ensembles
      sorted_pre = sort(pre_dists[i].egg_histogram,cols = :Age)
      sorted_post = sort(post_dists[i].egg_histogram,cols = :Age)
      pre = pre + (1.0*convert(Array,sorted_pre[:x1][1:25]) / pre_dists[i].n_samples)/ensembles
      post = post + (1.0*convert(Array,sorted_post[:x1][1:25]) / post_dists[i].n_samples)/ensembles
    end

#    pre =  convert(Array,pre_dists[1].egg_histogram[:x1][1:25]) / pre_dists[1].n_samples
#    post = convert(Array,post_dists[1].egg_histogram[:x1][1:25]) / post_dists[1].n_samples#
#    println()
#    println(sort(pre_dists[1].egg_histogram,cols = :Age))
#    println(pre)
#    println()
    #out = Output(statistics,pre,post)
    out = Output(statistics,pre,post)

    #return convert(Array,post_dists[1].age_histogram[:,2])
    return out
  end

  function run_mda(env::Environment,mda_pro::MDA_Program,p::MyParameters,dists::Epi_distributions)
    years = mda_pro.years
    sample_rate = years
    T = years*365
    ensembles = 250
    n_samples = div(T,sample_rate)

    v_pro = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")  # This is a "null" vaccination program, included because the functions above require one.
                                                                    # Would be better to have vaccination as an optional argument in funcs above

    statistics = initialise_statistics(n_samples,ensembles)

    #env = load_environment("environment_HIGHPREV1000.jld")
    #p.contact_rate = set_beta(mda_pro.R0,p,env,dists)
    statistics, pre_dists,post_dists = run_timesteps_MDA!(ensembles,T,sample_rate,env,statistics,p,mda_pro,v_pro)

    pre = zeros(25)
    post = zeros(25)

    for i in 1:ensembles
      sorted_pre = sort(pre_dists[i].egg_histogram,cols = :Age)
      sorted_post = sort(post_dists[i].egg_histogram,cols = :Age)
      pre = pre + (1.0*convert(Array,sorted_pre[:x1][1:25]) / pre_dists[i].n_samples)/ensembles
      post = post + (1.0*convert(Array,sorted_post[:x1][1:25]) / post_dists[i].n_samples)/ensembles
    end

    out = Output(statistics,pre,post)

    return out
  end

  # these next 2 functions are to set up parallel vaccination runs
  function vaccine_setup(dict::Dict{String,Any},p::MyParameters)
    programs = Array{Vaccination_Program,1}(0)

    for c in dict["intensity"]
      for i in 1:length(dict["ages"])
          for eff in  dict["efficacy"]
            push!(programs,Vaccination_Program(dict["duration"],dict["ages"][i],dict["coverages"][i],16,eff[1],eff[2],0.0,eff[3],c))
          end
      end
    end

    return programs
  end

  function vacmda_setup(dict::Dict{String,Any},p::MyParameters)

    programs = Array{VACMDA_Program,1}(0)
    for mda_cov in dict["MDA_cov"]
      for c in dict["intensity"]
        for i in 1:length(dict["ages"])
            for eff in  dict["efficacy"]
              push!(programs,VACMDA_Program(dict["duration"],dict["ages"][i],dict["coverages"][i],16,eff[1],eff[2],0.0,eff[3],c,mda_cov))
            end
        end
      end
    end

    return programs
  end

  function create_VACMDA()

    vacmda2 = Dict(
      "intensity" => ["low","high","medium"],
      "R0" => [6.88e-6],
      "duration" => 2.5,
      "MDA_cov" => [0.4],
      "efficacy" =>[[1.0,1.0,0.0]],
      "ages" => [[5,7,9,11]],
      "coverages" => [[0.6,0.6,0.6,0.7]]
    )

    vacmda5 = Dict(
      "intensity" => ["low","high","medium"],
      "R0" => [6.88e-6],
      "duration" => 5.0,
      "MDA_cov" => [0.4],
      "efficacy" =>[[1.0,1.0,0.0]],
      "ages" => [[1,6,11],[5,10,15]],
      "coverages" => [[0.85,0.6,0.7],[0.6,0.7,0.45]]
    )

    vacmda10 = Dict(
    "intensity" => ["low","high","medium"],
      "R0" => [6.88e-6],
      "duration" => 10.0,
      "MDA_cov" => [0.4],
      "efficacy" =>[[1.0,1.0,0.0]],
      "ages" => [[1,11],[5,15]],
      "coverages" => [[0.85,0.70],[0.60,0.45]]
    )

    vacmda20 = Dict(
    "intensity" => ["low","high","medium"],
      "R0" => [6.88e-6],
      "duration" => 20.0,
      "MDA_cov" => [0.4],
      "efficacy" =>[[1.0,1.0,0.0]],
      "ages" => [[1],[5]],
      "coverages" => [[0.85],[0.6]]
    )

    cases = Dict(
      "vacmda5" => vacmda5,
      "vacmda10" => vacmda10,
      "vacmda20" => vacmda20,

    )

    cases = [vacmda2,vacmda5,vacmda10,vacmda20]

    return cases
  end



  function create_programs()
    vaccine2 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 2.5,
    "ages" => [[5,7,9,11]],
    "coverages" => [[0.6,0.6,0.6,0.7]],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" =>[true,false],
    "MDA0" => [true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])

    vaccine5 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 5,
    "ages" => [[1,6,11],[5,10,15]],
    "coverages" => [[0.85,0.6,0.7],[0.6,0.7,0.45]],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" => [true,false],
    "MDA0" => [true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])

    vaccine10 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 10,
    "ages" => [[1,11],[5,15]],
    "coverages" => [[0.85,0.70],[0.60,0.45]],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" =>[true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])

    vaccine20 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 20,
    "ages" => [[1],[5]],
    "coverages" => [[0.85],[0.60]],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" =>[true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])

    cases = Dict(
    "vaccine2" => vaccine2,
    "vaccine5" => vaccine5,
    "vaccine10" => vaccine10,
    "vaccine20" => vaccine20)

    cases = [vaccine2,vaccine5,vaccine10,vaccine20]

    return cases
  end

  function create_extra_programs()

    sac = [5,6,7,8,9,10,11,12,13,14]
    sac_cov = [0.6,0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7,0.7]

    adult = collect(15:100)
    adult_cov = ones(adult)*0.45

    community = vcat(vcat([0,1,2,3,4],sac),adult)
    community_cov = vcat(vcat([0.85,0.85,0.6,0.6,0.6],sac_cov),adult_cov)


    vaccine5 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 5,
    "ages" => [sac,community],
    "coverages" => [sac_cov,community_cov],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" => [true,false],
    "MDA0" => [true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])

    vaccine10 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 10,
    "ages" => [sac,community],
    "coverages" => [sac_cov,community_cov],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" =>[true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])

    vaccine20 = Dict(
    "intensity" => ["low","high","medium"],
    "duration" => 20,
    "ages" => [sac,community],
    "coverages" => [sac_cov,community_cov],
    "R0" => [9.64e-6,6.88e-6],
    "N" => [500,1000,10000],
    "MassVaccinate0" =>[true,false],
    "efficacy" => [[0.75,0.75,0.0],[1.0,1.0,0.0]])


    cases = [vaccine5,vaccine10,vaccine20]

    return cases
  end

  function test_run(from_saved::Bool,save_env::Bool, b::Float64,years::Int,k::Float64,lambda::Float64,no::Int64)
    p = MyParameters()
    p.death_human = load_age_death_rates()
    #p.age_contact = load_contact_data("data/IietuneContact.txt")
    p.age_contact = load_contact_data("data/ICL.txt")
    #p.age_contact = load_contact_data("data/test.txt")
    p.gamma_r = Distributions.Gamma(k,1.0/k)
    p.max_fecundity = lambda

    if from_saved == true
      env = load_environment("environment1000.jld")
    else
      env = create_population(p,no,1000,50)

      println("created population")
      #update_env!(env,p)
    end

    distributions = get_distributions(env)
    r0 = calculateR0(p,env,distributions)
    println("R0 of original pop is ", r0)

    sample_rate = 25
    #years = 250
    T = years*365
    ensembles = 1
    n_samples = div(T,sample_rate)
    statistics = initialise_statistics(n_samples,ensembles)

    p.contact_rate = set_beta(b,p,env,distributions)
    #p.contact_rate = 0.000258

    stats, ee, dists = run_to_equilibrium!(env,1000,1,save_env,p,years);


    #stats = run_mda(env::Environment,mda_pro::MDA_Program,p::MyParameters,dists::Epi_distributions)

    println("runequilibrium success")
    println("mean aggregation is ",mean(stats.aggregation[end-1000:end,1]))
    println("mean prevalence is ",mean(stats.prevalence[end-1000:end,1]))
    println("end prevalence is ",stats.prevalence[end,1])
    println("mean H-prevalence is ",mean(mean(stats.stratified_prevalence[:,3,:],2)[end-1000:end]))
    println("end H-prevalence is ", mean(stats.stratified_prevalence[:,3,:],2)[end])
    println("contact rate is ", p.contact_rate)



    make_plots(stats,T,sample_rate)

#    '''
#    e = load_environment("environment1000.jld")
#
#    cases = create_programs()

#    for case in cases
#      println("duration is ", case["duration"])
#    end
#    case = cases[1]


#    # v_pro = Vaccination_Program(0.8,20,ages,40,0.9,2.5)
#    # s = run_vaccination(e,v_pro,p,distributions)
#    #make_plots(s,365*20,sample_rate)
#    '''
    return ee, dists
  end


  function find_stable_prevalence(bs::Array{Float64,1},years::Int)
    ps = zeros(bs)
    cs = zeros(bs)
    as = zeros(bs)
    rs = zeros(bs)
    iter = 0
    for b in bs
      iter += 1
      p = MyParameters()
      p.death_human = load_age_death_rates()
      #p.age_contact = load_contact_data("data/IietuneContact.txt")
      p.age_contact = load_contact_data("data/Iietune2.txt")

      env = create_population(p,5000,1000,50)
      println("created population")

      distributions = get_distributions(env)
      r0 = calculateR0(p,env,distributions)
      println("R0 of original pop is ", r0)

      sample_rate = 25
      #years = 250
      T = years*365
      ensembles = 5
      n_samples = div(T,sample_rate)
      statistics = initialise_statistics(n_samples,ensembles)

      p.contact_rate = set_beta(b,p,env,distributions)
      r0 = calculateR0(p,env,distributions)
      println("new R0 is ", r0)

      local stats, dists, ee
      for jj in 1:4
        env0 = deepcopy(env)
        stats, ee, dists = run_to_equilibrium!(env0,1000,1,false,p,years)
        if stats.prevalence[end,1] > 0.01
          break
        end
      end

      #stats = run_mda(env::Environment,mda_pro::MDA_Program,p::MyParameters,dists::Epi_distributions)

      println("runequilibrium success")
      #println("mean aggregation is ",mean(stats.aggregation[end-1000:end,1]))
      #println("mean prevalence is ",mean(stats.prevalence[end-1000:end,1]))
      #println("end prevalence is ",stats.prevalence[end,1])
      #println("mean H-prevalence is ",mean(mean(stats.stratified_prevalence[:,3,:],2)[end-1000:end]))
      #println("end H-prevalence is ", mean(stats.stratified_prevalence[:,3,:],2)[end])
      #println("contact rate is ", p.contact_rate)
      ps[iter] = mean(stats.prevalence[end-6000:end,1])
      cs[iter] = p.contact_rate
      as[iter] = mean(stats.aggregation[end-6000:end,1])
      rs[iter] = calculateR0(p,ee,dists)
    end
    return ps,cs,as,rs
  end

  function MDA_test(intensity::String,c::Float64)
    p = MyParameters()
    p.death_human = load_age_death_rates()
    p.age_contact = load_contact_data("data/ICL.txt")

    #p.gamma_r = Distributions.Gamma(k,1.0/k)
    #p.max_fecundity = lambda

    sample_rate = 16
    years = 16
    T = years*365
    ensembles = 250
    n_samples = div(T,sample_rate)

    ages = collect(1:1:150)
    cov = 1.0*zeros(ages)
    for i in 1:length(ages)
      if ages[i] < 5
        cov[i] = 0.0
      elseif ages[i] < 15
        cov[i] = c
      else
        cov[i] = 0
      end
    end
    mda_pro = MDA_Program(cov,years,0.86,0.0,3)

    if intensity == "high"

      env = load_environment("ENV1000HIGH5.200e-4.jld")
      p.contact_rate = 5.00e-4
      p.gamma_r = Distributions.Gamma(0.24,1.0/0.24)
      p.max_fecundity = 0.14
      mda_pro = MDA_Program(cov,years,0.86,0.0,1)

    end

    if intensity == "medium"
      env = load_environment("ENV1000MED1.543e-4.jld")
      p.contact_rate = 1.543e-4
      p.gamma_r = Distributions.Gamma(0.24,1.0/0.24)
      p.max_fecundity = 0.14
      mda_pro = MDA_Program(cov,years,0.86,0.0,2)
    end

    if intensity == "low"

      env = load_environment("ENV1000LOW1.390e-4.jld")
      p.contact_rate = 1.390e-4
      p.gamma_r = Distributions.Gamma(0.04,1.0/0.04)
      p.max_fecundity = 0.14
      mda_pro = MDA_Program(cov,years,0.86,0.0,3)

    end

    dists = get_distributions(env)
    #p.contact_rate = set_beta(r0,p,env,dists)


    v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")
    for i in 1:365*1
        update_env!(env,p,v)
    end

    out = run_mda(env,mda_pro,p,dists)
    save("mda_"*intensity*"_"*string(100*c)*"_snc0p4.jld", "out", out)
    #make_prevalence_plot(stats,T,sample_rate)
    #make_high_prevalence_plot(stats,T,sample_rate)

    return out
  end

  function generate_age_dist(from_saved::Bool)

    p = MyParameters()
    p.death_human = load_age_death_rates()
    p.age_contact = load_contact_data("data/ICL.txt")

    if from_saved == true
      env = load_environment("environment1000.jld")
    else
      env = create_population(p,1000,1000,50)
      println("created population")
      #update_env!(env,p)
    end

    distributions = get_distributions_byyear(env)
    v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")

    years = 500
    T = years*365

    #stats = run_mda(env::Environment,mda_pro::MDA_Program,p::MyParameters,dists::Epi_distributions)
    for i in 1:T
        update_env!(env,p,v)

        if (i%20 == 0) & (i> 15*T/20)
          update_distributions_byyear!(distributions,env)
        end
      end
    return distributions
  end

  function reintroduction(M::Int64,no::Int64)
    p = MyParameters()
    p.death_human = load_age_death_rates()
    p.age_contact = load_contact_data("data/Iietune2.txt")
    p.contact_rate = 7.759e-6
    env = load_environment("environment_MEDHIGH_con7p759e-6.jld")

    sample_rate = 16
    years = 16
    T = years*365
    ensembles = 50
    n_samples = div(T,sample_rate)
    env = load_environment("environment_MEDHIGH_con7p759e-6.jld")

    # Remove all worms
    for h in env.humans
      h.male_worms = zeros(p.worm_stages)
      h.female_worms = zeros(p.worm_stages)
      #h.male_worms = round(M/p.worm_stages)*ones(p.worm_stages)
      #h.female_worms = round(M/p.worm_stages)*ones(p.worm_stages)
      h.eggs = 0
      h.larvae = []
    end
    # remove worms from the environment
    env.larvae = []
    env.infective_larvae = 0


    #randomly select individual and infect with m worms
    #h_i = rand(1:length(env.humans))
    perm = randperm(length(env.humans))
    for h_i in perm[1:no]
      env.humans[h_i].male_worms = round(M/p.worm_stages)*ones(p.worm_stages)
      env.humans[h_i].female_worms = round(M/p.worm_stages)*ones(p.worm_stages)
    end

    distributions = get_distributions(env)
    statistics = initialise_statistics(n_samples,ensembles)
    v = Vaccination_Program(20,[1],[0.75],16,0.0,0.0,0.0,0.0,"high")
    env, d1, d2 = run_timesteps!(ensembles,T,sample_rate,env,statistics,p,v)

    return statistics, env, d2

  end

  # vaccination covesages by age, (from Stefano)
  function coverages(ages)
    covs = zeros(ages)
    for i in 1:length(ages)
      if i < 2
        covs[i] = 0.85
      elseif i < 10
        covs[i] = 0.6
      elseif i < 15
        covs[i] = 0.7
      else
        covs[i] = 0.45
      end
    end
    return covs
  end

  function efficacy_programs()

    programs = Array{Vaccination_Program,1}(0)
    durations = linspace(1,20,20)
    efficacies = 0:0.02:0.6
    for i in durations
      ages = collect(5.0:i:15.0)
      for j in efficacies
          ef = sqrt(1.0-j)
          push!(programs,Vaccination_Program(i,ages,coverages(ages),16,ef,ef,0.0,0.0,"high"))
      end
    end

    return programs
  end


#============================================================================#
end
