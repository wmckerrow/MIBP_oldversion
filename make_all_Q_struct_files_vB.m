function make_all_Q_struct_files_vB(strings, entropy_cutoff)
% This function makes one "Qstruct" file per vertex. These files store the
% partition function and number of structures for that vertex, for both the
% unconstrained space and the space with the entropy constraints applied.
% Thus, each file stores 4 values.
% INPUT: strings, a cell array of strings, where each string is a
%          "tree_path" of a vertex in the tree. e.g. '' is the root node, '0' is the
%           root node's left child, etc. (see B_master_script for details)
%        entropy_cutoff, a num, the entropy cutoff value used to constrain
%          base pairs with low entropy

program_constants;
answer = zeros(2, 2);
for i = 1:length(strings)
    disp('string:');
    disp(strings{i});
    if isempty(strings{i})
        pairs = [];
    else
        mibps = dlmread([RNA_NAME '_' strings{i} '_mibps.txt']);
        pairs = mibps(strings{i} == '1',:);
    end
    [forced_pairs,forced_nonpairs] = B_entropy_constraints3(strings{i}, entropy_cutoff, pairs);
    constrained_con = strcat(RNA_NAME, '_', strings{i}, '_entropy.CON'); %entropy-constrained constraint file
    create_constraint_file([], [], forced_pairs,forced_nonpairs, constrained_con);
    [energy_cons,~] = constrained_energy_and_entropy(constrained_con,seqfile);%entropy-constrained Q
    con_Q = energy2partition(energy_cons);
    disp('constrained Q');
    disp(con_Q);
    unconstrained_pfs = strcat(RNA_NAME, '_', strings{i}, '.pfs');
    uncon_Q = pfs2Q(unconstrained_pfs, 'temp.txt'); %regular Q
    disp('unconstrained Q');
    disp(uncon_Q);
        
    constrained_structs = count_structures_dave(constrained_con, seqfile); 
    unconstrained_con = strcat(RNA_NAME, '_', strings{i}, '.CON'); % regular constraint file
    
    unconstrained_structs = count_structures_dave(unconstrained_con, seqfile);
    answer(1, 1) = uncon_Q;
    answer(1, 2) = con_Q;
    answer(2, 1) = unconstrained_structs;
    answer(2, 2) = constrained_structs;
    output_file = strcat(RNA_NAME, '_', strings{i}, '_QandStruct.txt');
    dlmwrite(output_file, answer);
    %answer is 2x2 matrix to which results are saved
end

end

