library( data.table )

# the time before which you take the parents 
# there is an issue with right data censoring 
# not an issue in COVID-19 where infectious periods are not long but a problem in HIV
t_first =1980

ids <- fread('OnlyIds.csv')
# transimssions is dataframe which requires 3 columns: time_infected, ID_source and ID_recipient
trans = as.data.table(ids)

# get the IDs of the parent (i.e. in first period)
id_inf = trans[TimeOfInfection <= t_first   & TimeOfInfection , .(IdInfector = IdInfected ) ]
count_inf  = id_inf[,.N]

# get the number infected by each infected (put back those who infected zero)
t_trans = trans[ id_inf, on = "IdInfector", nomatch = 0][ ,.(N = .N), by = "IdInfector"]
t_trans = t_trans[ id_inf, on = "IdInfector"][ , .(IdInfector, N = ifelse( is.na(N), 0 , N))]
t = t_trans[ , .(count = .N), by = "N"]