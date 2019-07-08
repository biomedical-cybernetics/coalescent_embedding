MATLAB implementation of the angular separation index (ASI):
a quantitative measure to evaluate the separation of groups
over the circle circumference (2D) or sphere surface (3D).

The folder contains the main function "compute_angular_separation.m"
and an usage example ("run_example.m", "example_data.mat").


### REFERENCE ###

A. Muscoloni and C. V. Cannistraci (2019), "Angular separability of data clusters or network communities
in geometrical space and its relevance to hyperbolic embedding", arXiv:1907.00025


### INPUT ###

coords - 2D case: Nx1 vector containing for each sample the angular coordinates
         3D case: Nx2 matrix containing for each sample the coordinates (azimut,elevation)
         angular and azimuth coordinates must be in [0,2pi]
         elevation coordinates must be in [-pi/2,pi/2]
         the code automatically selects the 2D or 3D case depending on the size of the "coords" variable

labels - N labels for the samples indicating the group membership (numeric vector or cell of strings)

show_plot - [optional] 1 or 0 to indicate whether the plot of the results has to be shown or not (default = 1)

rand_reps - [optional] repetitions for evaluating random coordinates (default = 1000)

rand_seed - [optional] nonnegative integer seed for random number generator (by default a seed is created based on the current time)

worst_comp - [optional] 1 or 0 to indicate if the worst case should be approximated computationally or theoretically (default = 1); note that for the 3D case only the value 1 is valid

(NB: optional inputs not given or empty assume the default value)


### OUTPUT ###

index - overall index in [0,1], a value 1 indicates that all the groups
        are perfectly separated over the circle circumference (2D) or sphere surface (3D),
        the more the groups are mixed the more the index tends to 0, representing a worst-case scenario.

group_index - vector containing an index in [0,1] for each group,
              to assess its separation with respect to the other groups

pvalue - empirical p-value computed comparing the observed index with a null distribution
         of indexes obtained from random permutations of the coordinates
	 
	 
### CONTACT ###

For any problem, please contact:
Alessandro Muscoloni: alessandro.muscoloni@gmail.com
Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com