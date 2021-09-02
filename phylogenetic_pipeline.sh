#!/bin/bash
#
#Script to run a pipeline that generates genomic sequences from a run of the ibm model
#
# July 2021, F. Di Lauro

set -e #Make sure the script stops if any errors occur

# Main directory for running this code
# using $PWD as second input argument will use the current working directory
#num_lines is the number of lines to consider within the parameter f ile
#run_number is the actual line to read
#num_runs, is how many runs starting from run_number you are interested in 
main_dir=$1
num_lines=$2
counterfactual=0
run_number=$2
num_runs=0

model_dir='.'
fixed_param_dir='PARAMS_COMMUNITY1_ACCEPTED'
CoaTran_dir='CoaTran'
phylo_model='./coatran_transtree'
#(cd $CoaTran_dir; make clean; make;)

Output_dir='Outputs'

mkdir -p $Output_dir

#1 compile models
#(cd $model_dir/src; make clean; make all)

#2 make one run
echo "$main_dir"

(cd $model_dir/src; ./popart-simul.exe "$main_dir/$fixed_param_dir" $num_lines $counterfactual $run_number $num_runs)
cp "$main_dir/$fixed_param_dir/Output/"phylo*Run${run_number}* $main_dir/$Output_dir

#3 generate .favites out of transmission files
args_t='--start_date 1970 --end_date 2018'
args_trans_0="--inputfilename_trans_p0 $main_dir/$Output_dir/phylogenetic_transmission_CL01_Za_B_V1.2_patch0_Rand10_Run${run_number}_PCseed0_0.csv"
args_trans_1="--inputfilename_trans_p1 $main_dir/$Output_dir/phylogenetic_transmission_CL03_Za_C_V1.2_patch1_Rand10_Run${run_number}_PCseed0_0.csv"
args_ind_0="--inputfilename_indiv_p0 $main_dir/$Output_dir/phylogenetic_individualdata_CL01_Za_B_V1.2_patch0_Rand10_Run${run_number}_PCseed0_0.csv"
args_ind_1="--inputfilename_indiv_p1 $main_dir/$Output_dir/phylogenetic_individualdata_CL03_Za_C_V1.2_patch1_Rand10_Run${run_number}_PCseed0_0.csv"
args_other="--sampled_individuals 350 --start_sampling 2013 --end_sampling 2018"

python3 transmission_to_favites.py -o $main_dir/$Output_dir/output_Run${run_number} $args_t $args_trans_0 $args_trans_1 $args_ind_0 $args_ind_1    


#4 generate .newick out of favites through CoaTran
out_trans_net="$main_dir/$Output_dir/output_Run${run_number}_transmission_network.tsv"
out_samp_times="$main_dir/$Output_dir/output_Run${run_number}_sample_times.tsv" 
out_phylo_newick="$main_dir/$Output_dir/phylogenetic_Run${run_number}.newick"
(cd $CoaTran_dir; $phylo_model "$out_trans_net" "$out_samp_times" > "$out_phylo_newick";)


phylo_filename="$main_dir/$Output_dir/statistics_Run${run_number}.csv"
#5 Compute statistics from the phylogenetic tree 

Rscript phylostatistics.R $out_phylo_newick $phylo_filename 

#6 ????

#7 PROFIT
