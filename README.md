# TwoSpeciesOccupancy
2-spp occupancy model work

Paul,

I just uploaded a bunch of data and code.

A lot of it regards to the single species occupancy models, that I used to filter the covariates to be passed on to the multi species models.

My approach was to:
1. Develop single species occupancy models and retain the covariates whose 95% credible interval did not excluded zero (i.e. had a significant effect on occupancy or on detection);

2. Develop a series of asymmetrical multi species models, where I assess the effect of the dominant species on its respective subordinates (e.g. lynx as dominant over badgers, foxes, wildcats, martens and genets; then badgers as dominants of foxes, wildcats, martens and genets, and so on). This would give a total of 7 models (one for each dominant species).

3. Calculate the species interaction factor (see manuscript that I sent by e-mail) to assess the level of spatial avoidance/association for each species pair (and, if possible, the analogous metric for detection, which would provide a measure of probability co-detectability)

4. Try to relate the Species interactions factor with each species body mass ratio and trophic niche overlap, to try to understand if the strength in spatial avoidance is mediated by species relative body size and shared prey.

What do you think of this approach?
Does it look sound to you?

Given the high number of parameters that we need to estimate, I decided to model occupancy with a stacked approach. That is, bind the camera trapping data from the non breeding and breeding seasons together, and included a dummy covariate ‘season’ to account for occupancy changes between seasons. 
Because I am not interested in colonisation/extinction processes, but rather in covariate associations (more getting rid of their effects), I think this approach could allow us to get more precise estimates of their effects. What do you think?

Thanks,

Pedro
