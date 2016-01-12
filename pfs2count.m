function Q = pfs2Q(pfsfile, energyfile)
% This program converts from a pfs file to a count of the number of
% structures in the ensemble. Be sure to change the DATAPATH and
% ENERGY_COUNT to the correct paths -- that is, to the version of
% RNAstructure with energies set to zero.
% INPUTS: pfsfile, a string, the name of the pfs file.
%         energyfile, a string, the name of the energy file created as an
%           intermediate step. Since this file is created by this function,
%           energyfile can be any name.
%OUTPUT:  Q, a num, the "partition function" (i.e., since the energies are
%           zero, a count of the number of structures) of the ensemble.

program_constants;
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH_COUNT ' ; ' ENERGY_COUNT ' -s ' pfsfile '> ' energyfile]);
fid = fopen(energyfile, 'r');
energy_line = fgetl(fid);
fclose(fid);
start = find(energy_line == ':') + 2;
finish = regexp(energy_line, 'kcal/mol');
energy = str2num(energy_line(start:finish-1));
Q = energy2partition(energy);

end
