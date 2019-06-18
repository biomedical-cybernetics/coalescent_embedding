function [x, ids] = create_matrix(id1, id2, weights)

% Authors:
% Alessandro Muscoloni, 2017

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

%%% INPUT %%%
% id1, id2 - vectors of numbers or cells of strings such that id1(i) links to id2(i);
%            links are considered undirected and self-loops are removed.
% weights - [optional] weights of the links.
%
% NB: if the weights are not provided, duplicated and bidirectional links are discarded;
%     if the weights are provided, the weights of duplicated and
%     bidirectional links are added up.
%
%%% OUTPUT %%%
% x - if there is one connected component, it is the adjacency matrix;
%     if there are more connected components, it is a cell of adjacency
%     matrices, sorted by decreasing size.
% ids - if there is one connected component, it is a vector that indicates
%       for each node its original id; if there are more connected components,
%       it is a cell of vectors, in the same order as the adjacency matrices.

% check input
narginchk(2,3);
if ~(isvector(id1) && isvector(id2) && length(id1)==length(id2))
    error('id1 and id2 must be vectors of the same length.')
end
if nargin==3
    if ~(isvector(weights) && length(weights)==length(id1) && isnumeric(weights) && all(weights>0))
        error('weights must be a vector of positive numbers with the same length as id1.')
    end
    weighted = 1;
else
    weighted = 0;
end

E = length(id1);
display(['Edge pairs in input: ' num2str(E)])

% unique list of original ids
ids_complete = unique(union(id1,id2));
N = length(ids_complete);
display(['Unique ids in input: ' num2str(N)])

% map original ids to node ids
display('Mapping original ids to node ids...')
node1 = zeros(E,1);
node2 = zeros(E,1);
if isnumeric(ids_complete)
    if any(isnan(ids_complete)) || any(isinf(ids_complete))
        error('Input contains NaN or Inf')
    end
    for i = 1:N
        node1(id1==ids_complete(i)) = i;
        node2(id2==ids_complete(i)) = i;
    end
else
    if any(strcmp(ids_complete,''))
        error('Input contains empty strings')
    end
    for i = 1:N
        node1(strcmp(id1,ids_complete{i})) = i;
        node2(strcmp(id2,ids_complete{i})) = i;
    end
end

if ~weighted
    temp = unique([node1 node2], 'rows');
    node1 = temp(:,1); node2 = temp(:,2); clear temp;
    display(['Duplicated edge pairs removed: ' num2str(E - length(node1))])
end

% create the adjacency matrix
display('Creating adjacency matrix...')
if weighted
    x_complete = sparse(node1,node2,weights,N,N);
    display(['Self-loops removed: ' num2str(sum(x_complete(speye(size(x_complete))==1)>0))])
    x_complete(speye(size(x_complete))==1) = 0;
    x_complete = x_complete + x_complete';
    E = sum(sum(x_complete>0))/2;
else
    x_complete = sparse(node1,node2,1,N,N);
    x_complete(speye(size(x_complete))==1) = 0;
    x_complete = max(x_complete,x_complete');
    E = sum(sum(x_complete))/2;
    display(['Self-loops and bidirectional edge pairs removed: ' num2str(length(node1) - E)])
end
display(['Number of nodes: ' num2str(N)])
display(['Number of edges: ' num2str(E)])

% create one adjacency matrix for each connected component
display('Computing connected components...')
[ncc, comp] = graphconncomp(x_complete, 'Directed', false);
display(['Number of connected components: ' num2str(ncc)])
if ncc==1
    x = x_complete;
    ids = ids_complete;
else
    comp_size = zeros(ncc,1);
    for i = 1:ncc
        comp_size(i) = sum(comp==i);
    end
    [~,idx] = sort(comp_size, 'descend');
    x = cell(ncc,1);
    ids = cell(ncc,1);
    for i = 1:ncc
        x{i} = x_complete(comp==idx(i),comp==idx(i));
        ids{i} = ids_complete(comp==idx(i));
        E = sum(sum(x{i}>0))/2;
        display(['Component ' num2str(i) ': N=' num2str(comp_size(idx(i))) ' E=' num2str(E)])
    end
end