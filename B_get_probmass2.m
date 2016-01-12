function p = B_get_probmass2(tree_path, energy)
%This function returns the probability mass of the cluster defined by
%tree_path
% INPUTS: tree_path, a string of 0's and 1's representing the identity of
%           the cluster
%         energy is the total energy
% OUTPUT: p, a double between 0 and 1, the probability mass of the cluster

program_constants;

%create energy txt file
pfsname = strcat(RNA_NAME, '_', tree_path, '.pfs');
energy_file = strcat(RNA_NAME, '_', tree_path, '_energy.txt');
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' ENSEMBLE_ENERGY_PROG ' -s ' pfsname '> ' energy_file]);

%parse energy txt file
fid = fopen(energy_file, 'r');
energy_line = fgetl(fid);
fclose(fid);
start = find(energy_line == ':') + 2; % the value of the energy falls between indices "start" and "finish"
finish = regexp(energy_line, 'kcal/mol');
energy_cluster = str2num(energy_line(start:finish-1));

p = energy2partition(-energy+energy_cluster);

end

