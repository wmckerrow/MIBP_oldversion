function Q = energy2partition(energy)
% Converts from "ensemble energy" (see
% http://rna.urmc.rochester.edu/Text/EnsembleEnergy.html) to partition
% function, Q

RT = 0.61633;
Q = exp(energy/(-1*RT));
end

