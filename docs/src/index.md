# Schistoxpkg.jl

A package to run an individual based mode of a schistosomiasis outbreak. Generally people uptake larvae based on a contact rate defined by their age, along with
some predisposition which is chosen from a gamma distribution with mean 1, but some specified level of variance.

All parameters are stored in the parameters.jl file in the src folder.

The model has a parameter which defines the time step that we take forward each time. Several functions are then called each
time step which simulate the run of the outbreak. This is repeated until we reach a specified number of steps, usually corresponding
to stepping forward a chosen number of years into the future.

The standard approach is given by the following set of processes which all have their own function to execute them.
Each time step we advance the time of the simulation by the length of the time step and also add this time step to the age of each individual.
There is a chosen period at which contact rates are updated for each individual, where we check if someone has aged into a different age bracket, resulting if their
level of contact has changed.

We then calculate the total number of worms within each individual and the number of pairs of worms a person has.
These numbers are used to calculate how many eggs someone will produce. The number of eggs is chosen from a negative binomial distribution with mean equal to the
number of worm pairs multiplied by the max fecundity parameter and then multiplied by an exponential function which calculates the density dependent reduction in eggs produced,
`λ wp exp(-wp z)`.
We then kill the worms within human hosts at a given rate, which is based on average worm lifespan. 

Eggs are then hatched into the environment, with egg release dependent on the age specific contact rate of each individual.
Humans are given an age of death when they are born, which is based on some chosen death rates for each age group. We check each time step if anyone has outlived their age of death and if they have, they are then removed from the population.
Cercariae from the environment are then uptaken to each surviving individual based on their predisposition and contact rate. These immediately become worms within the human host.

We then perform any interventions which are due to take place at this point in time after which we will cull the miracidia and cercariae in the environment by a chosen percentage. After this we will add births to the population which occur at some specified rate.


There are other versions of this basic approach, where we don't age the population or include births and deaths and also where the population is aged but every death is simply matched with a birth, resulting in the population being kept constant.

```@index
```

```@autodocs
Modules = [Schistoxpkg]
```
