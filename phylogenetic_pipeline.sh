#!/bin/bash
#
#Script to run a pipeline that generates genomic sequences from a run of the ibm model
#
# July 2021, F. Di Lauro

set -e #Make sure the script stops if any errors occur

# Main directory for running this code
# using $PWD as second input argument will use the current working directory

main_dir=$1


#1 compile models

#2 make one run

#3 generate .favites out of transmission files

#4 generate .newick out of favites through CoaTran

#5 generate sequences out of .newick through phangorn

#6 optional: render the .dot file in a pdf
