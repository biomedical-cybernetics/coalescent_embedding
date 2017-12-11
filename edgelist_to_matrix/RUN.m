% Authors:
% Alessandro Muscoloni, 2017

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

% read the edge list file
display('Reading edge list file...')
edges = dlmread('opsahl_11.txt');

% create the adjacency matrix
[x, ids] = create_matrix(edges(:,1), edges(:,2));

% if there are multiple connected components
% take only the largest one
if iscell(x)
    display('Only the largest connected component is used.')
    x = x{1};
    ids = ids{1};
end

% save adjacency matrix and node ids
save('opsahl_11.mat', 'x', 'ids')