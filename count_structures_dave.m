function num_structs = count_structures_dave(constraint_file, seqfile)
% This function counts the number of structures in an ensemble. It does it
% using RNAstructure, by setting the energies to zero. This way, the
% recursions yield the number of structures instead of the partition
% function. You will need to create and compile a separate version of RNAstructure with
% the energies set to 0. Make sure to set PARTITION_COUNT to the
% appropriate path -- it should be the partition function of your
% zero-energy version of RNAstructure.
% INPUTS: constraint_file, a string, the name of the ensemble's constraint
%           file. if constraint_file equals 'none', then there are no
%           constraints
%         seqfile, a string, the name/path to (if it's not in the working
%           directory) the sequence file
% OUTPUT: num_structs, a num, the number of structures

program_constants;

pfsfile = 'temp_counting.pfs';

if strcmp(constraint_file, 'none')
    system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH_COUNT ' ; ' PARTITION_COUNT ' ' seqfile ' ' pfsfile ' >/dev/null']);
else
    system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH_COUNT ' ; ' PARTITION_COUNT ' -c ' constraint_file ' ' seqfile ' ' pfsfile ' >/dev/null']);
end
num_structs = pfs2count(pfsfile, 'temp.txt');

end

