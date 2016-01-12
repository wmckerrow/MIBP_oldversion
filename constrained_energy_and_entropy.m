function [energy,entropy] = constrained_energy_and_entropy(confile,seqfile)
%This function returns the probability mass and total base pair entropy of
%the cluster defined by tree_path.
% INPUTS: tree_path, a string of 0's and 1's representing the identity of
%           the cluster
% OUTPUT: total energy for that node

program_constants;

%create energy txt file
pfsname = 'temp.pfs';
energy_file = 'temp_energy.txt';
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PARTITION_PROG ' -c ' confile ' ' seqfile ' ' pfsname ' >/dev/null']);
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' ENSEMBLE_ENERGY_PROG ' -s ' pfsname '> ' energy_file]);

%parse energy txt file
fid = fopen(energy_file, 'r');
energy_line = fgetl(fid);
fclose(fid);
start = find(energy_line == ':') + 2; % the value of the energy falls between indices "start" and "finish"
finish = regexp(energy_line, 'kcal/mol');
energy = str2num(energy_line(start:finish-1));

%create probs file
probs_file = 'temp_prob.txt';
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PROB_PROG ' -t ' pfsname ' ' probs_file ' >/dev/null']);

%calulcate entropy
bp_probs = read_bp_probs('temp_prob.txt');
bp_ents = real(shannon_entropy([bp_probs(:,3) 1-bp_probs(:,3)]'))';
entropy = sum(bp_ents);

end