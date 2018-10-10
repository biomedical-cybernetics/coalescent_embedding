MATLAB implementation of the angular separation score:
a quantitative measure to evaluate the separation of groups
over the circle circumference (2D) or sphere surface (3D).

The folder contains the main function "compute_angular_separation.m"
and an usage example ("run_example.m", "example_data.mat").


### REFERENCE ###

Cacciola et al. (2017), "Coalescent embedding in the hyperbolic space
unsupervisedly discloses the hidden geometry of the brain", arXiv:1705.04192

NB: the mathematical formula adopted is the one named as "score_w2" in the Suppl. Algorithm 1 of the reference

### INPUT ###

coords - 2D case: Nx1 vector containing for each sample the angular coordinates
         3D case: Nx2 matrix containing for each sample the coordinates (azimut,elevation)
         angular and azimuth coordinates must be in [0,2pi]
         elevation coordinates must be in [-pi/2,pi/2]
         the code automatically selects the 2D or 3D case depending on the size of the "coords" variable

labels - N labels for the samples indicating the group membership (numeric vector or cell of strings)

show_plot - [optional] 1 or 0 to indicate whether the plot of the results has to be shown or not (default = 1)

rand_reps - [optional] repetitions for evaluating random coordinates (default = 1000)

rand_seed - [optional] seed for random number generator (default = 1)

(NB: optional inputs not given or empty assume the default value)


### OUTPUT ###

score - overall score in [0,1], a value 1 indicates that all the groups
        are perfectly separated over the circle circumference (2D) or sphere surface (3D),
        the more the groups are mixed the more the score tends to 0, representing a worst-case scenario.

group_scores - vector containing a score in [0,1] for each group,
               to assess its separation with respect to the other groups

pvalue - empirical p-value computed comparing the observed score with a null distribution
         of scores obtained from random permutations of the coordinates
	 
	 
### CONTACT ###

For any problem, please contact:
Alessandro Muscoloni: alessandro.muscoloni@gmail.com
Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com