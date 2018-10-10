% Usage example of the function "compute_angular_separation.m"

% load coordinates (azimuth,elevation) and the group labels for the samples
load('example_data.mat', 'coords_2D', 'coords_3D', 'labels')

% score 2D: evaluation of the separation over the circle circumference
[score_2D, group_scores_2D, pvalue_2D] = compute_angular_separation(coords_2D, labels);

% score 3D: evaluation of the separation over the sphere surface
[score_3D, group_scores_3D, pvalue_3D] = compute_angular_separation(coords_3D, labels);