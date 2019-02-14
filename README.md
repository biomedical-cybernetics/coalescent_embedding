# Coalescent embedding

## Reference

A. Muscoloni, J. M. Thomas, S. Ciucci, G. Bianconi, and C. V. Cannistraci, "Machine learning meets complex networks via coalescent embedding in the hyperbolic space", Nature Communications, 8:1615, 2017.
https://www.nature.com/articles/s41467-017-01825-5

## Folders description

* *coemb_svds_eigs*  
  It contains the main function for performing the coalescent embedding: *coalescent_embedding.m*.
  The dimension reduction techniques exploits the MATLAB function *eigs* and the function *lansvd* from the PROPACK library, which have a time complexity of O(kN^2). Since in our case k=2 or k=3, it is practically O(N^2).
  The other files are support functions for computing the edge betweenness centrality (from the MatlabBGL library) or for computing the SVD (from the PROPACK library).
  The subfolder *mex_all_platforms* contains the MEX files for different platforms.

  The support functions of the MatlabBGL library have been downloaded at:  
  http://mathworks.com/matlabcentral/fileexchange/10922-matlabbgl  
  The support functions of the PROPACK library have been downloaded at:  
  https://github.com/mavenlin/PropackMatlab4Windows  
  https://github.com/epfl-lts2/unlocbox/tree/master/test_bench/private

* *coemb_svd_eig*  
  It contains the main function for performing the coalescent embedding: *coalescent_embedding.m*.
  The dimension reduction techniques exploits the MATLAB built-in functions *svd* and *eig*, which have a time complexity of O(N^3).
  The other files are support functions for computing the edge betweenness centrality (from the MatlabBGL library).
  The subfolder *mex_all_platforms* contains the MEX files for different platforms.

* *visualization_and_evaluation*  
  - *plot_embedding.m*: plot the embedding given the network and the coordinates.
  - *angular_alignment.m*: align two embeddings of the same network, useful for a comparison of the plots.
  - *greedy_routing.m*: evaluate the greedy routing given the network and the coordinates.
  - *compare_embedding.m*: evaluate the comparison of two embeddings computing HD-correlation and C-score.
  - *angular_separation_index*: evaluate the separation of groups over the circle circumference or sphere surface (for details please see the README within the subfolder).

* *usage_example*  
  It contains a script *RUN_EXAMPLE.m* with two usage examples, all the functions needed for the embedding, visualization and evaluation (described at the previous points) and two example networks (*opsahl_11.mat* and *PSO_network.mat*).
  1. Example 1
     - 2D embedding of PSO network using RA1-LE-EA
     - evaluate comparison of embedding with original coordinates
     - plot comparison of embedding with original coordinates
  2. Example 2
     - 3D embedding of opsahl-11 network using RA1-ISO
     - evaluation of 3D greedy routing
     - 3D plot colored by ground-truth communities

* *edgelist_to_matrix*  
  It contains the function *create_matrix.m* for creating an adjacency matrix starting from an edge list, an example of edge list (*opsahl_11.txt*) and a script *RUN.m* with an usage example of the function.
	 
## Contact

For any problem, please contact:
* Alessandro Muscoloni: alessandro.muscoloni@gmail.com
* Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com