function Q = pfs2Q(pfsfile, energyfile)
% This function converts from a pfs file to the partition function value.
% See http://rna.urmc.rochester.edu/Text/partition.html for information on
% the pfs file.
% INPUTS: pfsfile, a string, the name of the pfs file
%         energyfile, a string, the name of the energy file that will be
%           saved as an intermediate step. Since this file is created
%           during the running of this function, energyfile can be any
%           name.
% OUTPUT: Q, a num, the partition function value.

program_constants;
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' ENSEMBLE_ENERGY_PROG ' -s ' pfsfile '> ' energyfile]);
fid = fopen(energyfile, 'r');
energy_line = fgetl(fid);
fclose(fid);
start = find(energy_line == ':') + 2;
finish = regexp(energy_line, 'kcal/mol');
energy = str2num(energy_line(start:finish-1));
Q = energy2partition(energy);

end
