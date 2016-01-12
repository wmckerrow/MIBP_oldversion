% This script is run as a first step in almost all the functions. It
% establishes variables such as the name and length of the RNA sequence,
% any bases or base pairs that are forced to be paired or unpaired, and the
% name of needed paths. Feel free to add things to this script as needed.

%**************************************************************************
%************************** MODEL PARAMETERS ******************************
%**************************************************************************

max_conf_splits = evalin('base','max_conf_splits'); % max number of conflicting basepairs to look for
hibp_cutoff = evalin('base','hibp_cutoff'); % probability need to be a highly probable base pair

stem_thresh = evalin('base','stem_thresh'); % By how many factors we allow the proabilities of bps in a stem to differ
min_cluster = evalin('base','min_cluster'); % Throw away clusters smaller than this
MI_cutoff = evalin('base','MI_cutoff'); % Stop when weighted MI falls below this cutoff

num_MIBPs = evalin('base','num_MIBPs'); % max number of MIBPs, a.k.a. the number of leaves in the tree - 1 
entropy_cutoff = evalin('base','entropy_cutoff'); % entropy cutoff used in "entropy constraints"

eta = evalin('base','eta');  % throw away pairs with prob less than eta or greater than 1-eta
             % when calculating MIBP (see luan lin's thesis for why) 
delta = evalin('base','delta');  % minimum bp prob

FORCED_PAIRS = evalin('base','FORCED_PAIRS'); 
FORCED_NONPAIRS = evalin('base','FORCED_NONPAIRS');
FORCED_UNPAIRED = evalin('base','FORCED_UNPAIRED');


%**************************************************************************
%*********************** PATHS ********************************************
%**************************************************************************

% path to specific RNAstructure programs
CT2DOT_PROG = evalin('base','CT2DOT_PROG');
DOT2CT_PROG = evalin('base','DOT2CT_PROG');
DATA_PATH = evalin('base','DATA_PATH');
PARTITION_PROG = evalin('base','PARTITION_PROG');
PROB_PROG = evalin('base','PROB_PROG');
STOCHASTIC_PROG = evalin('base','STOCHASTIC_PROG');
ENSEMBLE_ENERGY_PROG = evalin('base','ENSEMBLE_ENERGY_PROG');
DRAW_PROG = evalin('base','DRAW_PROG');
EFN_PROG = evalin('base','EFN_PROG');
FOLD_PROG = evalin('base','FOLD_PROG');
DATAPATH = evalin('base','DATAPATH');

% path to data tables for counting
PARTITION_COUNT = evalin('base','PARTITION_COUNT');
ENERGY_COUNT = evalin('base','ENERGY_COUNT');
DATAPATH_COUNT = evalin('base','DATAPATH_COUNT');

% path to gcc compiler
GCC = evalin('base','GCC');

% ********** NUMBER OF SAMPLES IN STOCHASTIC SAMPLING *********************

NUM_SAMPLES = evalin('base','NUM_SAMPLES');

%**************************************************************************
%******************************* RNA SEQUENCE *****************************
%**************************************************************************
RNA_NAME = evalin('base','RNA_NAME'); % What files created will be named
RNA_LENGTH = evalin('base','RNA_LENGTH'); % Length of the RNA molecule
seqfile = evalin('base','seqfile'); % Location of the seqeunce file (.fa or .seq)
native_structure = evalin('base','native_structure'); % Location of the native structure file (leave empty if unknown)
