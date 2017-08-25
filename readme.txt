(1) To run the executable, make this folder you current matlab directory and
run the master_script_hibp.m script. It runs modified RNAstructure executable,
which were copmiled for Linux (CentOS release 6.5). Using standard
RNAstructure executables is okay, but the probability of entropy constraints
will not be correct.

(2) To change which RNA molecule is being considered edit the program_constants.txt file.

(3) After completion, visualization.html will appear. Open this with firefox to view a visulization similar to Figure 3 of the paper. Click on a node of the tree to get a detailed view. Files for the visualization are located in the jsons folder. To retain a visualization simply rename visualization.html and make sure that it is in the same folder as a folder called 'jsons' containing the appropriate files.

(4) The results folder contains a large dump of data related to the run. If you are not interested in this, you may wish to delete it before running again. The format for the file names is:
<RNA_NAME>_<node_name>_<kind>
where kind is
energy: the total free energy of structures in the node (region or well in the paper). This value is used to calculate the probability of a region/well.
entropy: a constraint file used by RNAstructure to specificy the entropy constraints
probs: the probability plot output from RNAstructure containing log10 probability for each base pair with positive probability.
mibps: the base pairs whose presence/absence defines that node
Finally there is a <RNA_NAME>.mat file. If loaded in matlab it will give access to the summary statistics given in the paper:
RNA_LENGTH: length of the sequence
n: number of regions
probs: probility of regions without entropy constraints
probs_le: probility of regions with entropy constraints
all_structs: number of possible structure
structs: structures in each region without entropy constraints
structs_le: structures in each region with entropy constraints
sum_of_bp_entropy: sum of base pair entropy for the ensemble
entropies_le: sum of base pair entropy for each region with entropy constraints
centroids: centroids of each region
native_cluster: regions that contains the native structure
tree_strings: binary string name of node (1=go right, 0=go left)
