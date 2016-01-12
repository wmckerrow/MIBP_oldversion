function [mi] = mutual_info(p1, p2)
% calculate MI of two binary random variables from samples
% INPUTS: An array that gives the value of the first random variable in
% each sample, an array that gives the value of the second random variable in
% each sample.

program_constants;

% Count how many time each r.v. is 1
p1_count = sum(p1);
p2_count = sum(p2);

% Marginal for first r.v.
prob1 = [p1_count (NUM_SAMPLES - p1_count)]/NUM_SAMPLES;

% Marginal for second r.v.
prob2 = [p2_count (NUM_SAMPLES - p2_count)]/NUM_SAMPLES;

% Joint distribution
prob3 = [sum(p1==0 & p2==0) sum(p1==0 & p2==1) ...
         sum(p1==1 & p2==0) sum(p1==1 & p2==1)]/NUM_SAMPLES;
  
% I(X;Y) = H(X) + H(Y) - H(X,Y)
mi = shannon_entropy(prob1) + shannon_entropy(prob2) - ...
    shannon_entropy(prob3);

end