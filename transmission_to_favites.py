"""
Generate FAVITES files of transmission events from IBM output.  

This creates two dot files, 
1) who-infected-whom and when, with same name 
of output arg + _transmission_network.tsv

2) randomly individuals sampled 1 year after infection, 
same name of output file with suffix _sample_times 


Usage: 

python python/generate_dot.py \
    --start_date 1970 --end_date 1990\
    --outfilename ./output/transmissions_1970_1990\
    --inputfilename_trans_p0=$trans_patch0\
    --inputfilename_indiv_p0=$indiv_patch0\
    --inputfilename_trans_p1=$trans_patch1\
    --inputfilename_indiv_p1=$indiv_patch1\
    --sampled_individuals=50


runfile('transmission_to_favites.py', args = '-o output --start_date 1970 --end_date 2018 --inputfilename_trans_p0 phylogenetic_transmission_CL01_Za_B_V1.2_patch0_Rand10_Run1_PCseed0_0.csv --inputfilename_trans_p1 phylogenetic_transmission_CL03_Za_C_V1.2_patch1_Rand10_Run1_PCseed0_0.csv --inputfilename_indiv_p0 phylogenetic_individualdata_CL01_Za_B_V1.2_patch0_Rand10_Run1_PCseed0_0.csv --inputfilename_indiv_p1 phylogenetic_individualdata_CL03_Za_C_V1.2_patch1_Rand10_Run1_PCseed0_0.csv')    

"""
import pandas as pd, numpy as np
import argparse

# Separator in transmission files
SEP = ","
TAB = "\t"

# Red/blues
col_risk_low = "#2c7bb6" # (dark blue); " #abd9e9" (light blue)
col_risk_med = "#fdae61"
col_risk_high = "#d7191c"

# Blue/green for male/female (see popart_style_guide)
col_male = "#1f78b4"
col_female = "#e41a1c" #"#33a02c"

def consolidate_transmission(all_infection_events):
    """
    Combine two transmission dataframes into one file, ensuring that the list is 
    ordered in time, and that if A infects B and B infects C in the same time-slot,
    then theo order is A -> B, B->C and not B->C, A->B

    Parameters
    ----------
    dftrans_0, dftrans_1 : pandas.DataFrame
        Dataframes output by the IBM,  transmission files for patch 0 and 1 respectively


    Returns
    -------
    dftransmissions: pandas.DataFrame
        Dataframe of the same form of input, sorted by time

    """
    dftransmissions=all_infection_events.sort_values('TimeOfInfection',ignore_index=True)    
    #I could do groupby but I think it might mess up the order when going back?
    uniquetimes = dftransmissions.TimeOfInfection.unique()
    
    for time in uniquetimes:
        rows = dftransmissions.loc[dftransmissions.TimeOfInfection==time]
        infectors = rows.IDINFECTOR
        infected  = rows.IDINFECTED
        intersection = np.intersect1d(infectors,infected)  
        firstrow = rows.index[0]
        
        for element in intersection:    
            #Move intersection elements (when infected) to the top of the list

            #place_infector = dftransmissions.loc[dftransmissions.IDINFECTOR == element]
            place_infected = dftransmissions.loc[dftransmissions.IDINFECTED == element]
            temp_0 = dftransmissions.iloc[firstrow]
            dftransmissions.iloc[firstrow] = dftransmissions.iloc[place_infected.index[0]]
            dftransmissions.iloc[place_infected.index[0]] = temp_0
            firstrow = firstrow + 1
    
    dftransmissions = dftransmissions.reset_index(drop=True)

    #Need to repeat because I am not sure about the indices swapping
    for time in uniquetimes:
        rows = dftransmissions.loc[dftransmissions.TimeOfInfection==time]
        infectors = rows.IDINFECTOR
        infected  = rows.IDINFECTED
        intersection = np.intersect1d(infectors,infected)  
        firstrow = rows.index[0]
        
        for element in intersection:    
            place_infected = dftransmissions.loc[dftransmissions.IDINFECTED == element]
            infector= dftransmissions.iloc[place_infected.index[0]].IDINFECTOR      
            place_infector = dftransmissions.loc[dftransmissions.IDINFECTED == infector]
            
            if infector in intersection:
                if place_infected.index[0] < place_infector.index[0]:
                    #This means that the infector of the element has been infected before infecting
                    
                    #Therefore, swap them!
                    temp_0 = dftransmissions.iloc[firstrow]
                    dftransmissions.iloc[firstrow] = dftransmissions.iloc[place_infector.index[0]]
                    dftransmissions.iloc[place_infector.index[0]] = temp_0
                    firstrow = firstrow + 1                    
    
    
    return dftransmissions.reset_index(drop=True)   
    

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
    hj
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

      
    parser.add_argument('-n','--sampled_individuals',  required = False, \
        help = "number of sampled individuals", type = int, default = 50)

    parser.add_argument('-s','--start_sampling',  required = False, \
        help = "when to start sampling individuals", type = int, default = 1990)
    parser.add_argument('-e','--end_sampling',  required = False, \
        help = "when to end sampling individuals", type = int, default = 2018)
                
        
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
    dftransmissions = pd.concat([df_trans_patch0, df_trans_patch1],ignore_index=True)
    # Subset transmission events to only those of interest
    condition1 = dftransmissions.TimeOfInfection > args.start_date
    condition2 = dftransmissions.TimeOfInfection <= args.end_date
    #condition3 = dftransmissions.PATCHINFECTED.isin(args.patch)
    all_infection_events = dftransmissions[condition1 & condition2 ]
    
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

    all_infection_events = consolidate_transmission(all_infection_events)
    
    # Create the dot file with time as the x-axis
    f = open(args.outfilename + '_transmission_network.tsv', 'w')   
    g = open(args.outfilename + '_sample_times.tsv', 'w')

    #Sort by date

    
    f.write("None"+TAB+"SUPERFAKER"+TAB+"1968"+'\n')
    #g.write("SUPERFAKER"+TAB+"1970"+'\n')
    
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
    
      # If the infected person is a seed case, don't plot the transmission event.  
        if (row.IDINFECTOR != "-1_-1") and (row.IDINFECTOR != "-1_2"):
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


            #EDGE<TAB>u<TAB>v<TAB>attributes (csv or .)<TAB>(d)irected or (u)ndirected
                    
            f.write(row.IDINFECTOR+TAB+row.IDINFECTED+TAB+str(row.TimeOfInfection)+'\n')
            #g.write(row.IDINFECTED+TAB+str(row.TimeOfInfection+1)+'\n')
        else:
            fakecase = "Fake_%d"%index
            f.write("SUPERFAKER"+TAB+fakecase+TAB+"1969"+'\n')
            #g.write(fakecase+TAB+str(row.TimeOfInfection+1)+'\n')
           
            f.write(fakecase+TAB+row.IDINFECTED+TAB+str(row.TimeOfInfection)+'\n')
            #g.write(row.IDINFECTED+TAB+str(row.TimeOfInfection+1)+'\n')
                        
    f.close()
    

    condition1 = all_infection_events.TimeOfInfection >= args.start_sampling
    condition2 = all_infection_events.TimeOfInfection <= args.end_sampling
    
    
    
    tosample =  all_infection_events[condition1 & condition2]

    sampled_index = np.sort(np.random.choice(tosample.index, args.sampled_individuals, replace=False))
    
    sampled_infected = tosample.loc[sampled_index]
    for index, row in sampled_infected.iterrows():
        g.write(row.IDINFECTED+TAB+str(row.TimeOfInfection+1)+'\n')
    g.close()
    
    