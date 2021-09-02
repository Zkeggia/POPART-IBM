#!/usr/bin/env python3
"""
Script to merge transmission and individual files so that the complete transmission network for a
single transmission event is generated.  ID and patch are unique identifiers of individuals in the 
model so a new ID column is made for the individual file which consists of patch-specific ID and 
patch number: <Id>_<Patch>.  
Usage
python generate_complete_transmission_network.py \
    --patch0_dir
    --patch1_dir
    --community
    --arm
    --version
    --patch
    --Rand
    --Run
    --PCseed
    --CF
    --output_dir
Returns
-------
Two output files are saved to the <output_dir> input argument.  These filenames have the same prefix as the input files, so `phylogenetic_individualdata` and `phylogenetic_transmission` for individual and transmission files respectively.  All other arguments are kept the same except patch is changed to `patchall`.  
Filenames are therefore of the following form: 
    
    phylogenetic_individualdata_CL<community>_Za_<arm>_V<version>_patchall_Rand<Rand>_Run<Run>_PCseed<PCseed>_<CF_text>.csv
    phylogenetic_transmission_CL<community>_Za_<arm>_V<version>_patchall_Rand<Rand>_Run<Run>_PCseed<PCseed>_<CF_text>.csv
CF_text is determined from the "CF" argument.  Either "0" or "0_CF".  
W. Probert, 2019
"""

import sys, os, argparse
from os.path import join
import numpy as np, pandas as pd

trans_prefix="phylogenetic_transmission"
indiv_prefix="phylogenetic_individualdata"

if __name__ == "__main__":
    
    # Process the input argument
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--patch0_dir", type = str, help = "Directory of files for patch 0")
    
    parser.add_argument("--patch1_dir", type = str, help = "Directory of files for patch 1")
    
    parser.add_argument("--community", type = str, help = "Community number")
    
    parser.add_argument("--arm", type = str, help = "Arm of the community")
    
    parser.add_argument("--version", type = str, help = "IBM code version number", default = "1.2")
    
    parser.add_argument("--patch", type = str, help = "Patch number to treat as the inside patch (patch 0)", default = "0")
    
    parser.add_argument("--Rand", type = str, help = "Seed used", default = "1")
    
    parser.add_argument("--Run", type = str, help = "Run number", default = "0")
    
    parser.add_argument("--PCseed", type = str, 
        help = "Random seed used for the PC sampling process", default = "0")
    
    parser.add_argument("--CF", type = int, help = "Counterfactual or not (1 or 0 respectively)", default = 0)
    
    parser.add_argument("--output_dir", type = str, 
        help = "Directory to which CSV files should be written", default = "output")
    
    args = parser.parse_args()
    
    # Process input arguments for making file names
    community = args.community.zfill(2)
    
    if int(community) <= 12:
        country_text = "Za"
    else:
        country_text = "SA"
    
    # Adjust the counterfactual so that text is used instead of True/False
    CF_text = "0_CF" if args.CF == 1 else "0"

    # Adjust arm based upon CF data
    arm = "C" if args.CF == 1 else args.arm

    trans_file = "_".join([trans_prefix, "CL" + community, country_text, arm, \
        "V" + args.version, "patch" + args.patch, "Rand" + args.Rand, "Run" + args.Run, \
        "PCseed" + args.PCseed, CF_text]) + ".csv"
    
    indiv_file = "_".join([indiv_prefix, "CL" + community, country_text, arm, \
        "V" + args.version, "patch" + args.patch, "Rand" + args.Rand, "Run" + args.Run, \
        "PCseed" + args.PCseed, CF_text]) + ".csv"
    
    trans0 = pd.read_csv(join(args.patch0_dir, trans_file))
    trans1 = pd.read_csv(join(args.patch1_dir, trans_file))
    indiv0 = pd.read_csv(join(args.patch0_dir, indiv_file))
    indiv1 = pd.read_csv(join(args.patch1_dir, indiv_file))
    
    # Merge data on individuals from both patches into a full dataset
    # Assign a global ID to the data on individuals in each patch
    indiv0['ID'] = indiv0.Id.map(str) + "_0"
    indiv1['ID'] = indiv1.Id.map(str) + "_1"
    indiv0['PATCH'] = 0
    indiv1['PATCH'] = 1
    full_indiv = pd.concat([indiv0, indiv1])
    full_indiv['SEX'] = full_indiv['Sex']
    
    # Assign a global ID to the infected individuals in patch 0 and 1
    trans0['PATCH_INFECTED'] = 0
    trans1['PATCH_INFECTED'] = 1
    trans0['ID_INFECTED'] = trans0.IdInfected.map(str) + "_" + trans0.PATCH_INFECTED.map(str)
    trans1['ID_INFECTED'] = trans1.IdInfected.map(str) + "_" + trans1.PATCH_INFECTED.map(str)
    
    trans0['PATCH_INFECTOR'] = trans0['IsInfectorOutsidePatch']
    trans1['PATCH_INFECTOR'] = 1 - trans1['IsInfectorOutsidePatch']
    
    # Adjust PATCH_INFECTOR to record seed cases (otherwise they'll be set to have a patch of -2)
    trans1.loc[trans1.IsInfectorOutsidePatch == -1, 'PATCH_INFECTOR'] = -1

    trans0['ID_INFECTOR'] = trans0.IdInfector.map(str) + "_" + trans0.PATCH_INFECTOR.map(str)
    trans1['ID_INFECTOR'] = trans1.IdInfector.map(str) + "_" + trans1.PATCH_INFECTOR.map(str)
    full_trans = pd.concat([trans0, trans1])
    
    full_trans = pd.merge(full_trans, full_indiv[['ID', 'DoB', 'SEX']], left_on = 'ID_INFECTED', 
        right_on  ='ID', how = 'left', copy = False)
    
    full_trans = pd.merge(full_trans, full_indiv[['ID', 'DoB', 'SEX']], left_on = 'ID_INFECTOR', 
        right_on  ='ID', how = 'left')
    
    # Adjust the names of the columns in the transmission data set
    cols = full_trans.columns.values
    cols = ['DOB_INFECTED' if c == 'DoB_x' else c for c in cols]
    cols = ['DOB_INFECTOR' if c == 'DoB_y' else c for c in cols]
    cols = ['SEX_INFECTED' if c == 'SEX_x' else c for c in cols]
    cols = ['SEX_INFECTOR' if c == 'SEX_y' else c for c in cols]
    full_trans.columns = cols
    
    full_trans['AGE_INFECTED'] = full_trans.TimeOfInfection - full_trans.DOB_INFECTED
    full_trans['AGE_INFECTOR'] = full_trans.TimeOfInfection - full_trans.DOB_INFECTOR
    
    cols2save_trans = ['PATCH_INFECTED', 'PATCH_INFECTOR', 'ID_INFECTED', 'ID_INFECTOR', \
           'DOB_INFECTED', 'DOB_INFECTOR', 'TimeOfInfection', \
           'AGE_INFECTED', 'AGE_INFECTOR', 'SEX_INFECTED', 'SEX_INFECTOR', \
           'IdInfector', 'IdInfected',  'IsInfectorAcute', \
           'PartnerARTStatus', 'IsInfectorOutsidePatch', 'InfectorCD4', \
           'InfectorSPVL', 'InfectedSPVL', 'Infector_NPartners']
    
    cols2save_indiv = ["Id", "ID", "PATCH", "SEX", "DoB", "DoD", "HIV_pos", "RiskGp", \
        "t_diagnosed", "cd4_diagnosis", "cd4atfirstART", "t_1stARTstart", "t_1stVLsupp_start", \
        "t_1stVLsupp_stop"]
    
    # Calculate the number infected at the end of the period in question
    # in each age/sex category
    full_indiv = pd.merge(left = full_indiv, \
        right = full_trans[['ID_INFECTED', 'TimeOfInfection']], 
        left_on = 'ID', right_on = 'ID_INFECTED', how = 'left')
    
    # Pull out the filename stub for each of the files.  
    output_filename = "_".join(["CL" + community, country_text, \
        arm, "V" + args.version, "patchall", "Rand" + args.Rand, "Run" + args.Run, \
        "PCseed" + args.PCseed, CF_text]) + ".csv"
    
    full_indiv[cols2save_indiv].to_csv(join(args.output_dir, indiv_prefix + "_" + output_filename),\
        sep = ",", index = False)
    
    full_trans[cols2save_trans].to_csv(join(args.output_dir, trans_prefix + "_" + output_filename),\
        sep = ",", index = False)