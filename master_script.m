% This is the main script that does everything involved in dividing the
% posterior space of structures into a subset of mutually exclusive
% "clusters" based on the presence or absence of a subset of base pairs.

% I've been thinking of this (and describing it in the code) in terms of a
% rooted binary tree. The root vertex is the entire space of structures;
% its right child is the space containing one of the most informative base
% pairs (MIBPs), its left child is the space not containing that base pair. Each
% split of a node is based on a most informative base pair; the right child
% has that base pair, the left child doesn't.

%% Calculate partition function (Q) of the entire ensemble:
% (Q will be used to calculate the probability mass of the "sub-spaces")

pg_eval; % set program constants, such as RNA_NAME, FORCED_PAIRS, FORCED_NONPAIRS, etc.

system(['rm ' RNA_NAME '*MI.txt']);

confile = strcat(RNA_NAME, '_.CON'); % name of constraint file
create_constraint_file([], [], FORCED_PAIRS, FORCED_NONPAIRS, confile);
pfsfile = strcat(RNA_NAME, '_.pfs');
ctfile = strcat(RNA_NAME, '_.ct');

disp('sampling at root')
%system(['cp -f ' seqfile ' ~/RNA/RNAstructure/RNAstructure/exe/']);
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PARTITION_PROG ' -c ' confile ' ' seqfile ' ' pfsfile ' >/dev/null']);
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ;' STOCHASTIC_PROG ' ' pfsfile ' ' ctfile ' >/dev/null']);
% we include "export LD_LIBRARY_PATH..." because for some reason MATLAB
% sometimes uses the wrong C++ compiler

disp('calculating energy')
ensemble_energy_file = strcat(RNA_NAME, '__energy.txt');
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' ENSEMBLE_ENERGY_PROG ' -s ' pfsfile ' > ' ensemble_energy_file]);
fid = fopen(ensemble_energy_file, 'r');
energy_line = fgetl(fid); % this is a text line returned by RNAstructure. it contains the energy of the ensemble.
fclose(fid);
start = find(energy_line == ':') + 2; % start and finish find where in engergy_line the relevant number is.
finish = regexp(energy_line, 'kcal/mol');
energy = str2num(energy_line(start:finish-1));
Q = energy2partition(energy);

%% Find MIBP of entire ensemble (i.e., of the root vertex):
% disp('finding ensemble mibp')
% [mi pair] = B_get_MI('');
% MIfile = strcat(RNA_NAME, '__MI.txt');
% dlmwrite(MIfile, [mi pair]); %dlmwrite/dlmread used to read and write matlab matrices

%% Make the binary tree:
[tree_strings,all_tree_strings] = make_tree({''},{''}, num_MIBPs,MI_cutoff, energy);

% tree_strings is a cell array of strings, one string per leaf of the
% binary tree. A tree string of 001 corresponds to the leaf you get by
% starting at the root and walking one step left (0), one step left (0),
% and one step right (1). For instance, if the tree looks like this:

%            V1
%          /    \
%        V2     V3
%       /      /  \
%      V4     V5   V6

% then the leaves are V4, V5, and V6, and the tree strings are '00', '10',
% and '11', respectively.

% %% make all Q/struct files
% % These files store the number of structures and the partition functions of
% % each vertex in the tree. One file per vertex.
% 
% make_all_Q_struct_files_vB(all_tree_strings, entropy_cutoff);
% 
% %% make all json files
% make_all_json_vB(all_tree_strings);
% 
% %% move json files
% 
% system(['cp *' RNA_NAME '*json ../json_from_oscar/PDB_00426_2bit/']);
% system(['mv -f *' RNA_NAME '*json jsons/']);
% 
% %% save tree_strings and all_tree_strings
% tree_strings_file = strcat('B_',RNA_NAME,'_',num2str(num_MIBPs),'_tree_strings.txt');
% write_1D_cell(tree_strings,tree_strings_file);
% all_tree_strings_file = strcat('B_',RNA_NAME,'_',num2str(num_MIBPs),'_all_tree_strings.txt');
% write_1D_cell(all_tree_strings,all_tree_strings_file);

%% calculate total basepair entropy
make_all_probs_files_vB({''})
sum_of_bp_entropy = total_basepair_entropy('');

%% Calculate probablilities, structure and entropies
pair_probs_ents;

%% Find cluster of native structure
if ~isempty(native_structure)
    native_cluster = which_cluster3(native_structure,all_tree_strings);
end

%% Save everything
save([RNA_NAME '.mat']);

%% Calculations for visualization
make_and_store_jsons

%% Calculate conflicting bps
prob_of_conflicting

disp(['native_clust_prob:',num2str(probs(find(ismember(tree_strings,native_cluster))))])
disp(['exp_clust_prob:',num2str(sum(probs.^2))])
disp(['rand_clust_prob:',num2str(mean(probs))])
disp(['num_structs:',num2str(sum(structs)),',',num2str(sum(structs_le))])
disp(['prob_le:',num2str(sum(probs_le))])
disp(['entropy:',num2str(sum_of_bp_entropy),',',num2str(sum(probs.*entropies))])


%% Clean up
system('rm temp*');
