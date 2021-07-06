"""
Generate dot files of transmission events from IBM output.  

This creates two dot files, 1) that clusters singletons together, 2) that aligns nodes according 
to the time of infection (with the suffix "_time_aligned" on the end of the output file).  


Usage: 

python python/generate_dot.py \
    --start_date 1970 --end_date 1990\
    --outfilename ./output/transmissions_1970_1990\
    --inputfilename_trans_p0=$trans_patch0\
    --inputfilename_indiv_p0=$indiv_patch0\
    --inputfilename_trans_p1=$trans_patch1\
    --inputfilename_indiv_p1=$indiv_patch1\
    --fig_size 11.7 8.267\ 
    --colour_var=sex

The output files can then be processed:
> fdp -Tpdf output/transmissions_1970_1990.dot > graphics/transmissions_1970_1990.pdf


And the "time aligned" directed graphs need to be processed with the 'dot' engine:

> dot -Tpdf output/transmissions_1970_1990_time_aligned.dot > \
>     graphics/transmissions_1970_1990_time_aligned.pdf
"""

import pandas as pd, numpy as np
import argparse

# Separator in transmission files
SEP = ","
TAB = "    "

# Red/blues
col_risk_low = "#2c7bb6" # (dark blue); " #abd9e9" (light blue)
col_risk_med = "#fdae61"
col_risk_high = "#d7191c"

# Blue/green for male/female (see popart_style_guide)
col_male = "#1f78b4"
col_female = "#e41a1c" #"#33a02c"

def consolidate_individual_files(individual0, individual1):
    """
    Combine two individual dataframes into one file, creating a new "ID" column for 
    each individual by combining "Id" and "patch" columns with underscore separator.  
    
    Parameters
    ----------
    
    individual0, individual1: pandas.DataFrame
        Dataframes output by the IBM, individual files for patch 0 and 1 respectively.  
    
    Returns
    -------
    
    dfindividual: pandas.DataFrame
        Dataframe with new "ID" column that is the original "Id" and "patch" columns combined
    """
    
    individual0['patch'] = 0
    individual1['patch'] = 1
    
    # Create overarching individual file
    dfindividual = pd.concat([individual0, individual1])
    
    # Make a new ID column that combines patch-id and patch number.  
    dfindividual['ID'] = dfindividual['Id'].astype(str) + '_' + dfindividual['patch'].astype(str)
    
    return dfindividual


if __name__ == "__main__":
    # Process the input argument
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-n", "--graph_name", type = str, required = False, \
       help = "Graph name", default = "transmissions")
      
    parser.add_argument('-f','--fig_size', nargs = '+', required = False, \
        type = float, default = [19.2,10.8])
    
    parser.add_argument('-p','--patch', nargs = '+', required = False, \
        help = "Patch of infected individuals to focus on", type = int, default = [0])
    
    parser.add_argument('-c','--colour_var', required = False, \
        help = "Which variable to use for colouring nodes (risk or sex)", type = str, 
        default = "risk")
    
    parser.add_argument("--start_date", type = float, 
       help = "Start date at which to look at transmissions", required = True)
    
    parser.add_argument("--end_date", type = float, 
       help = "End date at which to look at transmissions", required = True)
    
    parser.add_argument("-o", "--outfilename", type = str, required = False, \
       help = "Output filename (excluding the filetype suffix)", default = 'output')
    
    parser.add_argument("--inputfilename_trans_p0", type = str, required = True, \
      help = "Input filename (transmission in patch 0)")

    parser.add_argument("--inputfilename_indiv_p0", type = str, required = True, \
      help = "Input filename (individuals in patch 0)")

    parser.add_argument("--inputfilename_trans_p1", type = str, required = True, \
      help = "Input filename (transmission in patch 1)")

    parser.add_argument("--inputfilename_indiv_p1", type = str, required = True, \
      help = "Input filename (individuals in patch 1)")
    
    args = parser.parse_args()
    
    df_indiv_patch0 = pd.read_csv(args.inputfilename_indiv_p0, sep = SEP)
    df_indiv_patch1 = pd.read_csv(args.inputfilename_indiv_p1, sep = SEP)
    df_trans_patch0 = pd.read_csv(args.inputfilename_trans_p0, sep = SEP)
    df_trans_patch1 = pd.read_csv(args.inputfilename_trans_p1, sep = SEP)
    
    # Consolidate the individual-level data
    dfindividual = consolidate_individual_files(df_indiv_patch0, df_indiv_patch1)
    
    df_trans_patch0['IDINFECTION'] = np.arange(df_trans_patch0.shape[0]) + 1
    df_trans_patch1['IDINFECTION'] = np.arange(df_trans_patch0.IDINFECTION.max(), \
        df_trans_patch0.IDINFECTION.max() + df_trans_patch1.shape[0]) + 1
    
    # All transmission events
    df_trans_patch0['PATCHINFECTED'] = 0
    df_trans_patch0['IDINFECTED'] = df_trans_patch0['IdInfected'].astype(str) + '_' + \
        df_trans_patch0['PATCHINFECTED'].astype(str)
    df_trans_patch0['PATCHINFECTOR'] = df_trans_patch0['IsInfectorOutsidePatch']
    df_trans_patch0['IDINFECTOR'] = df_trans_patch0['IdInfector'].astype(str) + '_' + \
        df_trans_patch0['PATCHINFECTOR'].astype(str)
        
    df_trans_patch1['PATCHINFECTED'] = 1
    df_trans_patch1['IDINFECTED'] = df_trans_patch1['IdInfected'].astype(str) + '_' + \
        df_trans_patch1['PATCHINFECTED'].astype(str)
    df_trans_patch1['PATCHINFECTOR'] = 1 - df_trans_patch1['IsInfectorOutsidePatch']
    df_trans_patch1['IDINFECTOR'] = df_trans_patch1['IdInfector'].astype(str) + '_' + \
        df_trans_patch1['PATCHINFECTOR'].astype(str)
    
    # Combine transmission files into a single file
    dftransmissions = pd.concat([df_trans_patch0, df_trans_patch1])
    
    # Subset transmission events to only those of interest
    condition1 = dftransmissions.TimeOfInfection > args.start_date
    condition2 = dftransmissions.TimeOfInfection <= args.end_date
    condition3 = dftransmissions.PATCHINFECTED.isin(args.patch)
    all_infection_events = dftransmissions[condition1 & condition2 & condition3]
    
    # Round to 4 DP.  
    all_infection_events.TimeOfInfection = all_infection_events.TimeOfInfection.round(decimals = 4)
    
    # Name columns that are to be kept/merged from the individual file
    vars2keep = ['ID', 'Sex', 'DoB', 'DoD', 'HIV_pos', 'RiskGp']
    
    temp = pd.merge(all_infection_events, dfindividual[vars2keep], 
        left_on = 'IDINFECTED', right_on = 'ID', how = 'left')
    
    # Rename the columns as being from 'INFECTED' individuals
    all_infection_events = temp.rename(columns = {'Sex':'SEXINFECTED', 'DoD': 'DODINFECTED',
        'DoB':'DOBINFECTED', 'HIV_pos': 'HIVPOSINFECTED', 'RiskGp': 'RISKGPINFECTED'})
    
    # Import RISK GROUP, DOB, DOD, HIV_POS of infector person
    temp = pd.merge(all_infection_events, dfindividual[vars2keep], 
        left_on = 'IDINFECTOR', right_on = 'ID', how = 'left')
    
    # Rename the columns as being INFECTED covariates
    all_infection_events = temp.rename(columns = {'Sex':'SEXINFECTOR', 'DoD': 'DODINFECTOR',
        'DoB':'DOBINFECTOR', 'HIV_pos': 'HIVPOSINFECTOR', 'RiskGp': 'RISKGPINFECTOR'})
    
    # Convert risk group variables to categorical ordinal variables
    for var in ['RISKGPINFECTOR', 'RISKGPINFECTED']:
        all_infection_events[var] = all_infection_events[var].astype('category')

    
    
    # Create the dot file with time as the x-axis
    f = open(args.outfilename + '_time_aligned.dot', 'w')
    
    f.write('digraph ' + args.graph_name + ' {\n')
    
    # Global graph attributes
    f.write(TAB + 'rankdir=LR;\n')
    f.write(TAB + 'ratio="fill";\n')
    f.write(TAB + 'overlap=false;\n')
    f.write(TAB + 'margin=0;\n')
    f.write(TAB + 'size="' + str(args.fig_size[0]) + "," + str(args.fig_size[1]) + '";\n')
    
    # Global node attributes
    f.write(TAB + 'node [shape = circle, width=1.5, ')
    f.write('label="", fixedsize=true, overlap=scalexy, splines=true];\n')#, label=""
    
    # Global edge attributes
    f.write(TAB + 'edge [arrowhead = normal, label="", color="#919191"];\n')
    
    f.write(TAB + "{\n")
    f.write(TAB + "node [style=invis];\n")
    f.write(TAB + "edge [style=invis];\n")
    
    f.write(TAB + str(all_infection_events.TimeOfInfection.min()))
    
    for i, t in enumerate(np.sort(all_infection_events.TimeOfInfection.unique())[1:]):
        f.write(" -> " + str(t))
        if (i % 4) != 0:
            f.write("\n" + TAB)
    f.write(";")
    f.write(TAB + "}\n")
    
    for t in np.sort(all_infection_events.TimeOfInfection.unique()):
        
        sub = all_infection_events.loc[all_infection_events.TimeOfInfection == t]
        f.write('{ rank = same; "' + str(t) + '"; ')
        for index, row in sub.iterrows():
            f.write('"' + str(row.IDINFECTED) + '"; ')
        f.write('}\n')
    
    for index, row in all_infection_events.iterrows():
        
        if row.SEXINFECTED == row.SEXINFECTOR:
            print("Error, same sex infection")
            print("IDINFECTOR", row.IDINFECTOR)
            print("IDINFECTED", row.IDINFECTED)
            break
        
        if args.colour_var == "risk":
            if row.RISKGPINFECTED == "L":
                col = col_risk_low
            elif row.RISKGPINFECTED == "M":
                col = col_risk_med
            else: 
                col = col_risk_high
        elif args.colour_var == "sex":
            if row.SEXINFECTED == "M":
                col = col_male
            else:
                col = col_female
        else: 
            # If not being coloured by sex or risk group then use light slate grey for node colour
            col = "#778899"
        
        f.write(TAB+TAB + '"' + row.IDINFECTED + '"' + \
            '[style=filled, color=seashell4, fillcolor="' + col + '"];\n')
        
        # If the infected person is a seed case, don't plot the transmission event.  
        if (row.IDINFECTOR != "-1_-1"):
            if args.colour_var == "risk":
                if row.RISKGPINFECTOR == "L":
                    col = col_risk_low
                elif row.RISKGPINFECTOR == "M":
                    col = col_risk_med
                else: 
                    col = col_risk_high
            elif args.colour_var == "sex":
                if row.SEXINFECTOR == "M":
                    col = col_male
                else:
                    col = col_female
            else: 
                # If not being coloured by sex or risk group then use light slate grey for 
                # node colour
                col = "#778899"
            
            f.write(TAB+TAB + '"' + row.IDINFECTOR + '"' + '[style=filled, ');
            f.write(' fillcolor="' + col + '"');
            #if row.IsInfectorOutsidePatch == 1:
            #    f.write(', penwidth=8, color=black');
            #else:
            f.write(', color=seashell4');
            f.write('];\n')
            
            #for (infector, infected) in infections:
            f.write(TAB+TAB + '"' + row.IDINFECTOR + '"' + ' -> ' + '"' + row.IDINFECTED + '"')
        
            # If this is an infection that occurred during the acute phase, then make the pen width
            # darker than if it had occurred during the chronic phase.  
            if row.IsInfectorAcute == 1:
                f.write(' [penwidth=30, color="#3B3B3B"]')
            else:
                f.write(' [penwidth=20]')
            f.write(';\n')
    f.write('}\n')
    f.close()
    
    
    # Patch 0 infection events.  
    # Subset the transmission data to the date of interest.  
    condition1 = (df_trans_patch0.TimeOfInfection >= args.start_date)
    condition2 = (df_trans_patch0.TimeOfInfection <= args.end_date)
    infection_events = df_trans_patch0[condition1 & condition2]
    
    # Create an ID that includes the patch
    infection_events['PATCHINFECTED'] = 0
    infection_events['IDINFECTED'] = infection_events['IdInfected'].astype(str) + '_' + \
        infection_events['PATCHINFECTED'].astype(str)
    
    infection_events['PATCHINFECTOR'] = infection_events['IsInfectorOutsidePatch']
    infection_events['IDINFECTOR'] = infection_events['IdInfector'].astype(str) + '_' + \
        infection_events['PATCHINFECTOR'].astype(str)
    
    # Name columns that are to be kept/merged from the individual file
    vars2keep = ['ID', 'Sex', 'DoB', 'DoD', 'HIV_pos', 'RiskGp']
    
    temp = pd.merge(infection_events, dfindividual[vars2keep], 
        left_on = 'IDINFECTED', right_on = 'ID', how = 'left')
    
    # Rename the columns as being from 'INFECTED' individuals
    infection_events = temp.rename(columns = {'Sex':'SEXINFECTED', 'DoD': 'DODINFECTED',
        'DoB':'DOBINFECTED', 'HIV_pos': 'HIVPOSINFECTED', 'RiskGp': 'RISKGPINFECTED'})
    
    # Import RISK GROUP, DOB, DOD, HIV_POS of infector person
    temp = pd.merge(infection_events, dfindividual[vars2keep], 
        left_on = 'IDINFECTOR', right_on = 'ID', how = 'left')
    
    # Rename the columns as being INFECTED covariates
    infection_events = temp.rename(columns = {'Sex':'SEXINFECTOR', 'DoD': 'DODINFECTOR',
        'DoB':'DOBINFECTOR', 'HIV_pos': 'HIVPOSINFECTOR', 'RiskGp': 'RISKGPINFECTOR'})
    
    # Convert risk group variables to categorical ordinal variables
    for var in ['RISKGPINFECTOR', 'RISKGPINFECTED']:
        infection_events[var] = infection_events[var].astype('category')

    
    # Find singletons (not involved in any transmission events during the period in question)
    condition1 = (dfindividual.HIV_pos == 1)
    condition2 = (dfindividual.DoB < args.start_date)
    condition3 = (dfindividual.DoD > args.end_date)
    condition4 = (~dfindividual.ID.isin(infection_events.IDINFECTOR))
    condition5 = (~dfindividual.ID.isin(infection_events.IDINFECTED))
    condition6 = (dfindividual.patch == 0)
    singletons = dfindividual[condition1 & condition2 & condition3 & \
        condition4 & condition5 & condition6]
    
    # Create the dot file
    f = open(args.outfilename + '.dot', 'w')
    
    f.write('digraph ' + args.graph_name + ' {\n')
    #f.write(TAB + 'rankdir=LR;\n')
    f.write(TAB + 'size="' + str(args.fig_size[0]) + "," + str(args.fig_size[1]) + '";\n')
    
    # Global node attributes
    f.write(TAB + 'node [shape = circle, label="", overlap=scalexy, splines=true];\n')
    f.write(TAB + 'graph [fontsize = 144, resolution=300];\n')
    
    # Global edge attributes
    f.write(TAB + 'edge [arrowhead = normal, label="", color="#919191"];\n')
    
    f.write(TAB + 'subgraph cluster_transmissions {\n')
    f.write(TAB+TAB + 'label="Transmissions";\n')
    f.write(TAB+TAB + 'color=white;\n')
    for index, row in infection_events.iterrows():
        
        if row.SEXINFECTED == row.SEXINFECTOR:
            print("Error, same sex infection")
            print("IDINFECTOR", row.IDINFECTOR)
            print("IDINFECTED", row.IDINFECTED)
            break
        
        if row.RISKGPINFECTED == "L":
            col = col_risk_low
        elif row.RISKGPINFECTED == "M":
            col = col_risk_med
        else: 
            col = col_risk_high
        
        f.write(TAB+TAB + '"' + row.IDINFECTED + '"' + \
            '[style=filled, color=seashell4, fillcolor="' + col + '"];\n')
        
        if row.RISKGPINFECTOR == "L":
            col = col_risk_low
        elif row.RISKGPINFECTOR == "M":
            col = col_risk_med
        else: 
            col = col_risk_high
        
        f.write(TAB+TAB + '"' + row.IDINFECTOR + '"' + '[style=filled, ');
        f.write(' fillcolor="' + col + '"');
        #if row.IsInfectorOutsidePatch == 1:
        #    f.write(', penwidth=8, color=black');
        #else:
        f.write(', color=seashell4');
        
        f.write('];\n')
        
        #for (infector, infected) in infections:
        f.write(TAB+TAB + '"' + row.IDINFECTOR + '"' + ' -> ' + '"' + row.IDINFECTED + '"')
        
        # If this is an infection that occurred during the acute phase, then make the pen width
        # much thicker than a normal infection.  
        if row.IsInfectorAcute == 1:
            f.write(' [penwidth=8, color="#3B3B3B"]')
        else:
            f.write(' [penwidth=5]')
        
        f.write(';\n')
    
    f.write('}\n')
    
    f.write(TAB + 'subgraph cluster_singletons {\n')
    f.write(TAB+TAB + 'label="Singletons";\n')
    f.write(TAB+TAB + 'color=white;\n')
    for index, row in singletons.iterrows():
        if row.RiskGp == "L":
            col = col_risk_low
        elif row.RiskGp == "M":
            col = col_risk_med
        else: 
            col = col_risk_high
        
        f.write(TAB+TAB + '"' + row.ID + '"' + \
            '[style=filled, color=seashell4, fillcolor="' + col + '"];\n')
    
    f.write(TAB + '}\n')
    f.write('}\n')
    f.close()

