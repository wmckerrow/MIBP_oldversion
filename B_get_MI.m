function [max_mi max_pair] = B_get_MI(tree_path)
% This function finds the most informative base pair of a cluster.
% INPUT: tree_path, a string of 0's and 1's, the identity of the cluster
% OUTPUTS: max_mi, a double, the maximum average mutual information of a
%           base pair in the cluster
%          max_pair, the base pair having average mutual information (with
%            all other base pairs) equal to max_mi. I.e., this is the most
%            informative base pair

program_constants;

% sample structures
ctname = strcat(RNA_NAME, '_', tree_path, '.ct');
structures = readct(ctname, RNA_LENGTH);

% get mutual information
[max_mi max_pair avg_mi] = rna_mi3(structures);

end

