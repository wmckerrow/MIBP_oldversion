% Calculate useful infomation about the resulting regions: probability,
% forced pairs, entropy and number of strucutures with and without entropy
% constraints.

%% Without entropy constraints
program_constants;
disp('Calculating cluster info')
n = length(tree_strings)
probs = zeros(n,1);
pairs = {};
entropies = zeros(n,1);
structs = zeros(n,1);
for i=1:length(tree_strings)
    cluster = tree_strings{i};
    if length(cluster) > 0
        mibps = dlmread([RNA_NAME '_' cluster '_mibps.txt']);
        pairs{i} = mibps(cluster == '1',:); % A cell listing the pairs present in each cluster
    else
        pairs{i} = zeros(0,2)
    end
    probs(i) = B_get_probmass2(cluster,energy); % Probability of each cluster
    entropies(i) = real(total_basepair_entropy(cluster)); % Sum of base pair entropy for each cluster
    create_constraint_file([], [], pairs{i}, [], 'temp.CON'); 
    structs(i) = count_structures_dave('temp.CON', seqfile); % Number of structures in each cluster
end

%% With entropy constraints
program_constants;
n = length(tree_strings)
disp('calculating low entropy info')
probs_le = zeros(n,1);
entropies_le = zeros(n,1);
structs_le = zeros(n,1);
for i=1:length(tree_strings)
    i
    cluster = tree_strings{i};
    [forced_pairs,forced_nonpairs] = B_entropy_constraints3(cluster,entropy_cutoff,pairs{i});
    create_constraint_file([], [], forced_pairs,forced_nonpairs, 'temp.CON');
    [energy_cons,entropies_le(i)] = constrained_energy_and_entropy('temp.CON',seqfile);
    probs_le(i) = energy2partition(-energy+energy_cons);
    structs_le(i) = count_structures_dave('temp.CON', seqfile);
end

%% Count all structures in ensemble
create_constraint_file([], [], [],[], 'temp.CON');
all_structs = count_structures_dave('temp.CON', seqfile);
